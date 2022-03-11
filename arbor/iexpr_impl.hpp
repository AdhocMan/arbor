#pragma once

// Implementations for inhomogeneous expressions.

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

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

struct distance : public iexpr_interface {
    distance(double s, mlocation l) : scale(s), loc(std::move(l)) {}

    virtual double eval(const mprovider &p, const mcable &c) const override {
        auto eval_loc = mlocation{c.branch, (c.dist_pos + c.prox_pos) / 2};

        if(loc.branch == eval_loc.branch)
            return scale * std::abs(p.embedding().integrate_length(eval_loc, loc));

        // Locations on different branches.
        // Find first common parent branch. Branch id of parent is always smaller.
        msize_t eval_b = eval_loc.branch;
        msize_t loc_b = loc.branch;
        while(eval_b != loc_b) {
            if(eval_b > loc_b)
                eval_b = p.morphology().branch_parent(eval_b);
            else
                loc_b = p.morphology().branch_parent(loc_b);
        }

        // If mnpos, locations are on different sides of root. Take distance to root in this case.
        // Otherwise, take distance to end of parent branch
        const auto base_loc = eval_b == mnpos ? mlocation{0, 0.0} : mlocation{eval_b, 1.0};

        // compute distance to distal end of parent branch and add together
        return scale * (std::abs(p.embedding().integrate_length(loc, base_loc)) +
               std::abs(p.embedding().integrate_length(eval_loc, base_loc)));
    }

    double scale;
    mlocation loc;
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
