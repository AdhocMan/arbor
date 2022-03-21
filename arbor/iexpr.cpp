// Implementations for inhomogeneous expressions.

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <variant>

#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/arbexcept.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/iexpr.hpp>

namespace arb {

namespace iexpr_impl {
namespace {
struct scalar : public iexpr_interface {
    scalar(double v) : value(v) {}

    virtual double eval(const mprovider &, const mcable &) const override {
        return value;
    }

    double value;
};

struct radius : public iexpr_interface {
    radius(double s) : scale(s) {}

    virtual double eval(const mprovider &p, const mcable &c) const override {
        auto eval_loc = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};
        return scale * p.embedding().radius(eval_loc);
    }

    double scale;
};

template<class> inline constexpr bool always_false_v = false;

struct distance : public iexpr_interface {
    distance(double s, std::variant<mlocation_list, mextent> l) : scale(s), locations(std::move(l)) {}

    virtual double eval(const mprovider &p, const mcable &c) const override {
        auto eval_loc = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};

        return scale * std::visit([&](auto&& arg) -> double {
            using T = std::decay_t<decltype(arg)>;

            double min_dist = std::numeric_limits<double>::max();
            if constexpr (std::is_same_v<T, mlocation_list>) {
                for(const auto& loc : arg) {
                    min_dist = std::min(min_dist, compute_distance(p, eval_loc, loc));
                }
            } else if constexpr (std::is_same_v<T, mextent>) {
                for(const auto& c : arg) {
                    // distance is 0 if location within extent
                    if (c.branch == eval_loc.branch && c.prox_pos < eval_loc.pos && c.dist_pos > eval_loc.pos)
                        return 0.0;
                    min_dist = std::min(min_dist, compute_distance(p, eval_loc, prox_loc(c)));
                    min_dist = std::min(min_dist, compute_distance(p, eval_loc, dist_loc(c)));
                }
            } else {
                static_assert(always_false_v<T>, "non-exhaustive visitor!");
            }
            return min_dist;
            }, locations);
    }

    double compute_distance(const mprovider &p, const mlocation& loc_a, const mlocation& loc_b) const {

        if(loc_a.branch == loc_b.branch)
            return std::abs(p.embedding().integrate_length(loc_a, loc_b));

        // Locations on different branches.
        // Find first common parent branch. Branch id of parent is always smaller.
        msize_t branch_a = loc_a.branch;
        msize_t branch_b = loc_b.branch;
        while(branch_a != branch_b) {
            if(branch_b == mnpos || (branch_a != mnpos && branch_a > branch_b))
                branch_a = p.morphology().branch_parent(branch_a);
            else
                branch_b = p.morphology().branch_parent(branch_b);
        }

        // If mnpos, locations are on different sides of root. Take distance to root in this case.
        // Otherwise, take distance to end of parent branch
        const auto base_loc = branch_a == mnpos ? mlocation{0, 0.0} : mlocation{branch_a, 1.0};

        // compute distance to distal end of parent branch and add together
        return std::abs(p.embedding().integrate_length(loc_a, base_loc)) +
               std::abs(p.embedding().integrate_length(loc_b, base_loc));
    }

    double scale;
    std::variant<mlocation_list, mextent> locations;
};

struct add : public iexpr_interface {
    add(std::unique_ptr<iexpr_interface> l, std::unique_ptr<iexpr_interface> r) :
        left(std::move(l)), right(std::move(r)) {}

    virtual double eval(const mprovider &p, const mcable &c) const override {
        return left->eval(p, c) + right->eval(p, c);
    }

    std::unique_ptr<iexpr_interface> left;
    std::unique_ptr<iexpr_interface> right;
};

struct mul : public iexpr_interface {
    mul(std::unique_ptr<iexpr_interface> l, std::unique_ptr<iexpr_interface> r) :
        left(std::move(l)), right(std::move(r)) {}

    virtual double eval(const mprovider &p, const mcable &c) const override {
        return left->eval(p, c) * right->eval(p, c);
    }

    std::unique_ptr<iexpr_interface> left;
    std::unique_ptr<iexpr_interface> right;
};

} // namespace
} // namespace iexpr_impl


iexpr iexpr::scalar(double value) {
    return iexpr(iexpr_type::scalar, std::make_tuple(value));
}

iexpr iexpr::distance(double scale, locset loc) {
    return iexpr(iexpr_type::distance, std::make_tuple(scale, std::variant<locset, region>(std::move(loc))));
}

iexpr iexpr::distance(double scale, region reg) {
    return iexpr(iexpr_type::distance, std::make_tuple(scale, std::variant<locset, region>(std::move(reg))));
}

iexpr iexpr::radius(double scale) {
    return iexpr(iexpr_type::radius, std::make_tuple(scale));
}

iexpr iexpr::diameter(double scale) {
    return iexpr(iexpr_type::diameter, std::make_tuple(scale));
}

iexpr iexpr::add(iexpr left, iexpr right) {
    return iexpr(iexpr_type::add, std::make_tuple(std::move(left), std::move(right)));
}

iexpr iexpr::mul(iexpr left, iexpr right) {
    return iexpr(iexpr_type::mul, std::make_tuple(std::move(left), std::move(right)));
}


std::unique_ptr<iexpr_interface> thingify(const iexpr& expr, const mprovider& m) {
    switch (expr.type()) {
        case iexpr_type::scalar:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::scalar(
                        std::get<0>(std::any_cast<const std::tuple<double>&>(expr.args()))));
        case iexpr_type::distance: {
                const auto& scale = std::get<0>(std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));
                const auto& var = std::get<1>(std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));

                return std::visit([&](auto&& arg) {
                        return std::unique_ptr<iexpr_interface>(new iexpr_impl::distance(
                                scale, thingify(arg, m)));
                    }, var);
            }
        case iexpr_type::radius:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::radius(
                        std::get<0>(std::any_cast<const std::tuple<double>&>(expr.args()))));
        case iexpr_type::diameter:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::radius(
                        2.0 * std::get<0>(std::any_cast<const std::tuple<double>&>(expr.args()))));
        case iexpr_type::add:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::add(
                        thingify(std::get<0>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m),
                        thingify(std::get<1>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m)));
        case iexpr_type::mul:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::mul(
                        thingify(std::get<0>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m),
                        thingify(std::get<1>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m)));
    }
    return nullptr;
}


} // namespace arb
