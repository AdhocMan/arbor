#pragma once

// Implementations for inhomogeneous expressions.

#include <algorithm>
#include <cmath>
#include <limits>

#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/arbexcept.hpp>

namespace arb {
struct iexpr_interface {

    virtual double eval(const mprovider &p, const mcable &c) const = 0;

    virtual ~iexpr_interface() = default;
};

namespace iexpr_impl {
struct identity : public iexpr_interface {
    virtual double eval(const mprovider &, const mcable &) const override {
        return 1.0;
    }
};

struct distance : public iexpr_interface {
    distance(mlocation_list l) : locations(std::move(l)) {
        if (locations.empty())
            throw arbor_exception(
                "iexpr distance: Empty locset or region now allowed.");
    }

    virtual double eval(const mprovider &p, const mcable &c) const override {
        auto dist = std::numeric_limits<double>::max();
        std::for_each(locations.begin(), locations.end(), [&](const auto &loc) {
              dist = std::min(
                  std::abs(p.embedding().integrate_length(
                      loc, mlocation{c.branch, (c.dist_pos - c.prox_pos) / 2})), dist);
        });

        return dist;
    }

    mlocation_list locations;
  };
}

} // namespace arb
