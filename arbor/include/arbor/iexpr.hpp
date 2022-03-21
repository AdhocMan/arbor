#pragma once

// Implementations for inhomogeneous expressions.

#include <memory>
#include <any>

#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/primitives.hpp>

namespace arb {

enum class iexpr_type {
    scalar,
    distance,
    proximal_distance,
    distal_distance,
    radius,
    diameter,
    add,
    mul
};

struct iexpr {

    static iexpr scalar(double value);

    static iexpr distance(double scale, locset loc);

    static iexpr distance(double scale, region reg);

    static iexpr proximal_distance(double scale, locset loc);

    static iexpr proximal_distance(double scale, region reg);

    static iexpr distal_distance(double scale, locset loc);

    static iexpr distal_distance(double scale, region reg);

    static iexpr radius(double scale);

    static iexpr diameter(double scale);

    static iexpr add(iexpr left, iexpr right);

    static iexpr mul(iexpr left, iexpr right);

    iexpr_type type() const { return type_; }

    const std::any& args() const { return args_; }

private:
    iexpr(iexpr_type type, std::any args) : type_(type), args_(std::move(args)) {}

    iexpr_type type_;
    std::any args_;
};

struct iexpr_interface {

    virtual double eval(const mprovider &p, const mcable &c) const = 0;

    virtual ~iexpr_interface() = default;
};

std::unique_ptr<iexpr_interface> thingify(const iexpr& expr, const mprovider& m);

} // namespace arb
