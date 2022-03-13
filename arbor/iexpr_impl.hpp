#pragma once

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

namespace arb {

namespace iexpr_impl {
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

        return std::visit([&](auto&& arg) -> double {
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
            if(branch_a > branch_b)
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

} // namespace iexpr_impl

} // namespace arb
