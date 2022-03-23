// Implementations for inhomogeneous expressions.

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <variant>

#include <arbor/arbexcept.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/iexpr.hpp>
#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/util/any_visitor.hpp>

namespace arb {

namespace iexpr_impl {
namespace {

msize_t common_parent_branch(msize_t branch_a, msize_t branch_b, const morphology& m) {
    // Locations on different branches.
    // Find first common parent branch. Branch id of parent is
    // always smaller.
    while (branch_a != branch_b) {
        if (branch_b == mnpos || (branch_a != mnpos && branch_a > branch_b))
            branch_a = m.branch_parent(branch_a);
        else
            branch_b = m.branch_parent(branch_b);
    }

    return branch_a;
}

// Compute the distance in proximal direction. Will return nullopt if loc_prox is not between origin
// and loc_dist
std::optional<double> compute_proximal_distance(const mlocation& loc_prox,
    const mlocation& loc_dist,
    const mprovider& p) {

    // check order if on same branch
    if (loc_prox.branch == loc_dist.branch && loc_prox.pos > loc_dist.pos) return std::nullopt;

    // check if loc_prox branch is in proximal direction from loc_dist
    auto b = loc_dist.branch;
    while (b > loc_prox.branch) {
        b = p.morphology().branch_parent(b);
        if (b == mnpos) return std::nullopt;
    }
    if (b != loc_prox.branch) return std::nullopt;

    return p.embedding().integrate_length(loc_prox, loc_dist);
};

struct scalar: public iexpr_interface {
    scalar(double v): value(v) {}

    double eval(const mprovider&, const mcable&) const override { return value; }

    double value;
};

struct radius: public iexpr_interface {
    radius(double s): scale(s) {}

    double eval(const mprovider& p, const mcable& c) const override {
        auto loc_eval = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};
        return scale * p.embedding().radius(loc_eval);
    }

    double scale;
};

struct distance: public iexpr_interface {
    distance(double s, std::variant<mlocation_list, mextent> l):
        scale(s),
        locations(std::move(l)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        auto loc_eval = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};

        auto compute_distance =
            [](const mprovider& p, const mlocation& loc_a, const mlocation& loc_b) {
                if (loc_a.branch == loc_b.branch)
                    return std::abs(p.embedding().integrate_length(loc_a, loc_b));

                // If mnpos, locations are on different sides of root. Take
                // distance to root in this case. Otherwise, take distance to
                // end of parent branch
                const auto base_branch =
                    common_parent_branch(loc_a.branch, loc_b.branch, p.morphology());
                const auto base_loc =
                    base_branch == mnpos ? mlocation{0, 0.0} : mlocation{base_branch, 1.0};

                // compute distance to distal end of parent branch and add
                // together
                return std::abs(p.embedding().integrate_length(loc_a, base_loc)) +
                       std::abs(p.embedding().integrate_length(loc_b, base_loc));
            };

        return scale *
               std::visit(
                   arb::util::overload(
                       [&](const mlocation_list& arg) {
                           double min_dist = std::numeric_limits<double>::max();
                           for (const auto& loc: arg) {
                               min_dist = std::min(min_dist, compute_distance(p, loc_eval, loc));
                           }
                           return min_dist;
                       },
                       [&](const mextent& arg) {
                           double min_dist = std::numeric_limits<double>::max();
                           for (const auto& c: arg) {
                               // distance is 0 if location within extent
                               if (c.branch == loc_eval.branch && c.prox_pos < loc_eval.pos &&
                                   c.dist_pos > loc_eval.pos)
                                   return 0.0;
                               min_dist =
                                   std::min(min_dist, compute_distance(p, loc_eval, prox_loc(c)));
                               min_dist =
                                   std::min(min_dist, compute_distance(p, loc_eval, dist_loc(c)));
                           }
                           return min_dist;
                       }),
                   locations);
    }

    double scale;
    std::variant<mlocation_list, mextent> locations;
};

struct proximal_distance: public iexpr_interface {
    proximal_distance(double s, std::variant<mlocation_list, mextent> l):
        scale(s),
        locations(std::move(l)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        auto loc_eval = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};

        return scale *
               std::visit(arb::util::overload(
                              [&](const mlocation_list& arg) {
                                  std::optional<double> min_dist;
                                  for (const auto& loc: arg) {
                                      auto dist = compute_proximal_distance(loc_eval, loc, p);
                                      if (dist)
                                          min_dist = std::min(
                                              min_dist.value_or(std::numeric_limits<double>::max()),
                                              dist.value());
                                  }
                                  return min_dist.value_or(0.0);
                              },
                              [&](const mextent& arg) {
                                  std::optional<double> min_dist;
                                  for (const auto& c: arg) {
                                      if (c.branch == loc_eval.branch &&
                                          c.prox_pos < loc_eval.pos && c.dist_pos > loc_eval.pos)
                                          return 0.0;
                                      const mlocation loc{c.branch, c.dist_pos};
                                      auto dist = compute_proximal_distance(loc_eval, loc, p);
                                      if (dist)
                                          min_dist = std::min(
                                              min_dist.value_or(std::numeric_limits<double>::max()),
                                              dist.value());
                                  }
                                  return min_dist.value_or(0.0);
                              }),
                   locations);
    }

    double scale;
    std::variant<mlocation_list, mextent> locations;
};

struct distal_distance: public iexpr_interface {
    distal_distance(double s, std::variant<mlocation_list, mextent> l):
        scale(s),
        locations(std::move(l)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        auto loc_eval = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};

        return scale *
               std::visit(arb::util::overload(
                              [&](const mlocation_list& arg) {
                                  std::optional<double> min_dist;
                                  for (const auto& loc: arg) {
                                      auto dist = compute_proximal_distance(loc, loc_eval, p);
                                      if (dist)
                                          min_dist = std::min(
                                              min_dist.value_or(std::numeric_limits<double>::max()),
                                              dist.value());
                                  }
                                  return min_dist.value_or(0.0);
                              },
                              [&](const mextent& arg) {
                                  std::optional<double> min_dist;
                                  for (const auto& c: arg) {
                                      if (c.branch == loc_eval.branch &&
                                          c.prox_pos < loc_eval.pos && c.dist_pos > loc_eval.pos)
                                          return 0.0;
                                      const mlocation loc{c.branch, c.prox_pos};
                                      auto dist = compute_proximal_distance(loc, loc_eval, p);
                                      if (dist)
                                          min_dist = std::min(
                                              min_dist.value_or(std::numeric_limits<double>::max()),
                                              dist.value());
                                  }
                                  return min_dist.value_or(0.0);
                              }),
                   locations);
    }

    double scale;
    std::variant<mlocation_list, mextent> locations;
};

struct add: public iexpr_interface {
    add(std::unique_ptr<iexpr_interface> l, std::unique_ptr<iexpr_interface> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        return left->eval(p, c) + right->eval(p, c);
    }

    std::unique_ptr<iexpr_interface> left;
    std::unique_ptr<iexpr_interface> right;
};

struct sub: public iexpr_interface {
    sub(std::unique_ptr<iexpr_interface> l, std::unique_ptr<iexpr_interface> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        return left->eval(p, c) - right->eval(p, c);
    }

    std::unique_ptr<iexpr_interface> left;
    std::unique_ptr<iexpr_interface> right;
};

struct mul: public iexpr_interface {
    mul(std::unique_ptr<iexpr_interface> l, std::unique_ptr<iexpr_interface> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        return left->eval(p, c) * right->eval(p, c);
    }

    std::unique_ptr<iexpr_interface> left;
    std::unique_ptr<iexpr_interface> right;
};

struct div: public iexpr_interface {
    div(std::unique_ptr<iexpr_interface> l, std::unique_ptr<iexpr_interface> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        return left->eval(p, c) / right->eval(p, c);
    }

    std::unique_ptr<iexpr_interface> left;
    std::unique_ptr<iexpr_interface> right;
};

struct exp: public iexpr_interface {
    exp(std::unique_ptr<iexpr_interface> v): value(std::move(v)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        return std::exp(value->eval(p, c));
    }

    std::unique_ptr<iexpr_interface> value;
};

struct log: public iexpr_interface {
    log(std::unique_ptr<iexpr_interface> v): value(std::move(v)) {}

    double eval(const mprovider& p, const mcable& c) const override {
        return std::log(value->eval(p, c));
    }

    std::unique_ptr<iexpr_interface> value;
};

}  // namespace
}  // namespace iexpr_impl

iexpr::iexpr(double value) { *this = iexpr::scalar(value); }

iexpr iexpr::scalar(double value) { return iexpr(iexpr_type::scalar, std::make_tuple(value)); }

iexpr iexpr::distance(double scale, locset loc) {
    return iexpr(
        iexpr_type::distance, std::make_tuple(scale, std::variant<locset, region>(std::move(loc))));
}

iexpr iexpr::distance(double scale, region reg) {
    return iexpr(
        iexpr_type::distance, std::make_tuple(scale, std::variant<locset, region>(std::move(reg))));
}

iexpr iexpr::proximal_distance(double scale, locset loc) {
    return iexpr(iexpr_type::proximal_distance,
        std::make_tuple(scale, std::variant<locset, region>(std::move(loc))));
}

iexpr iexpr::proximal_distance(double scale, region reg) {
    return iexpr(iexpr_type::proximal_distance,
        std::make_tuple(scale, std::variant<locset, region>(std::move(reg))));
}

iexpr iexpr::distal_distance(double scale, locset loc) {
    return iexpr(iexpr_type::distal_distance,
        std::make_tuple(scale, std::variant<locset, region>(std::move(loc))));
}

iexpr iexpr::distal_distance(double scale, region reg) {
    return iexpr(iexpr_type::distal_distance,
        std::make_tuple(scale, std::variant<locset, region>(std::move(reg))));
}

iexpr iexpr::radius(double scale) { return iexpr(iexpr_type::radius, std::make_tuple(scale)); }

iexpr iexpr::diameter(double scale) { return iexpr(iexpr_type::diameter, std::make_tuple(scale)); }

iexpr iexpr::add(iexpr left, iexpr right) {
    return iexpr(iexpr_type::add, std::make_tuple(std::move(left), std::move(right)));
}

iexpr iexpr::sub(iexpr left, iexpr right) {
    return iexpr(iexpr_type::sub, std::make_tuple(std::move(left), std::move(right)));
}

iexpr iexpr::mul(iexpr left, iexpr right) {
    return iexpr(iexpr_type::mul, std::make_tuple(std::move(left), std::move(right)));
}

iexpr iexpr::div(iexpr left, iexpr right) {
    return iexpr(iexpr_type::div, std::make_tuple(std::move(left), std::move(right)));
}

iexpr iexpr::exp(iexpr value) { return iexpr(iexpr_type::exp, std::make_tuple(std::move(value))); }

iexpr iexpr::log(iexpr value) { return iexpr(iexpr_type::log, std::make_tuple(std::move(value))); }

std::unique_ptr<iexpr_interface> thingify(const iexpr& expr, const mprovider& m) {
    switch (expr.type()) {
        case iexpr_type::scalar:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::scalar(
                std::get<0>(std::any_cast<const std::tuple<double>&>(expr.args()))));
        case iexpr_type::distance: {
            const auto& scale = std::get<0>(
                std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));
            const auto& var = std::get<1>(
                std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));

            return std::visit(
                [&](auto&& arg) {
                    return std::unique_ptr<iexpr_interface>(
                        new iexpr_impl::distance(scale, thingify(arg, m)));
                },
                var);
        }
        case iexpr_type::proximal_distance: {
            const auto& scale = std::get<0>(
                std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));
            const auto& var = std::get<1>(
                std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));

            return std::visit(
                [&](auto&& arg) {
                    return std::unique_ptr<iexpr_interface>(
                        new iexpr_impl::proximal_distance(scale, thingify(arg, m)));
                },
                var);
        }
        case iexpr_type::distal_distance: {
            const auto& scale = std::get<0>(
                std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));
            const auto& var = std::get<1>(
                std::any_cast<const std::tuple<double, std::variant<locset, region>>&>(expr.args()));

            return std::visit(
                [&](auto&& arg) {
                    return std::unique_ptr<iexpr_interface>(
                        new iexpr_impl::distal_distance(scale, thingify(arg, m)));
                },
                var);
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
        case iexpr_type::sub:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::sub(
                thingify(std::get<0>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m),
                thingify(std::get<1>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m)));
        case iexpr_type::mul:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::mul(
                thingify(std::get<0>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m),
                thingify(std::get<1>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m)));
        case iexpr_type::div:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::div(
                thingify(std::get<0>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m),
                thingify(std::get<1>(std::any_cast<const std::tuple<iexpr, iexpr>&>(expr.args())), m)));
        case iexpr_type::exp:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::exp(
                thingify(std::get<0>(std::any_cast<const std::tuple<iexpr>&>(expr.args())), m)));
        case iexpr_type::log:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::log(
                thingify(std::get<0>(std::any_cast<const std::tuple<iexpr>&>(expr.args())), m)));
    }

    throw std::runtime_error("thingify iexpr: Unknown iexpr type");
    return nullptr;
}

}  // namespace arb
