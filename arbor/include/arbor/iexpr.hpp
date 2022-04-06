#pragma once

// Implementations for inhomogeneous expressions.

#include <any>
#include <memory>
#include <string>
#include <ostream>

#include <arbor/morph/primitives.hpp>
#include <arbor/morph/region.hpp>
#include <arbor/morph/locset.hpp>

namespace arb {

struct mprovider;

enum class iexpr_type {
    scalar,
    distance,
    proximal_distance,
    distal_distance,
    interpolation,
    radius,
    diameter,
    add,
    sub,
    mul,
    div,
    exp,
    log,
    named
};

struct iexpr {
    // Convert to scalar expr type
    iexpr(double value);

    iexpr_type type() const { return type_; }

    const std::any& args() const { return args_; }

    static iexpr scalar(double value);

    static iexpr distance(double scale, locset loc);

    static iexpr distance(double scale, region reg);

    static iexpr proximal_distance(double scale, locset loc);

    static iexpr proximal_distance(double scale, region reg);

    static iexpr distal_distance(double scale, locset loc);

    static iexpr distal_distance(double scale, region reg);

    static iexpr interpolation(double prox_value,
        locset prox_list,
        double dist_value,
        locset dist_list);

    static iexpr interpolation(double prox_value,
        region prox_list,
        double dist_value,
        region dist_list);

    static iexpr radius(double scale);

    static iexpr diameter(double scale);

    static iexpr add(iexpr left, iexpr right);

    static iexpr sub(iexpr left, iexpr right);

    static iexpr mul(iexpr left, iexpr right);

    static iexpr div(iexpr left, iexpr right);

    static iexpr exp(iexpr value);

    static iexpr log(iexpr value);

    static iexpr named(std::string name);


private:
    iexpr(iexpr_type type, std::any args): type_(type), args_(std::move(args)) {}

    iexpr_type type_;
    std::any args_;
};

std::ostream& operator<<(std::ostream& o, const iexpr& e);

inline iexpr operator+(iexpr a, iexpr b) { return iexpr::add(std::move(a), std::move(b)); }

inline iexpr operator-(iexpr a, iexpr b) { return iexpr::sub(std::move(a), std::move(b)); }

inline iexpr operator*(iexpr a, iexpr b) { return iexpr::mul(std::move(a), std::move(b)); }

inline iexpr operator/(iexpr a, iexpr b) { return iexpr::div(std::move(a), std::move(b)); }

inline iexpr operator+(iexpr a) { return a; }

inline iexpr operator-(iexpr a) { return iexpr::mul(-1.0, std::move(a)); }

struct iexpr_interface {

    virtual double eval(const mprovider& p, const mcable& c) const = 0;

    virtual ~iexpr_interface() = default;
};

std::shared_ptr<iexpr_interface> thingify(const iexpr& expr, const mprovider& m);

}  // namespace arb
