#pragma once

#include <arbor/cable_cell_param.hpp>
#include <arbor/common_types.hpp>
#include <arbor/export.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/util/lexcmp_def.hpp>

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace arb {

using network_hash_type = std::uint64_t;

struct ARB_SYMBOL_VISIBLE network_site_info {
    cell_gid_type gid;
    cell_kind kind;
    cell_tag_type label;
    mlocation location;
    mpoint global_location;

    ARB_ARBOR_API friend std::ostream& operator<<(std::ostream& os, const network_site_info& s);
};

ARB_DEFINE_LEXICOGRAPHIC_ORDERING(network_site_info,
    (a.gid, a.kind, a.label, a.location, a.global_location),
    (b.gid, a.kind, b.label, b.location, b.global_location))

struct ARB_SYMBOL_VISIBLE network_connection_info {
    network_site_info src, dest;
    double weight, delay;

    ARB_ARBOR_API friend std::ostream& operator<<(std::ostream& os, const network_connection_info& s);
};

ARB_DEFINE_LEXICOGRAPHIC_ORDERING(network_connection_info, (a.src, a.dest), (b.src, b.dest))

struct network_selection_impl;

struct network_value_impl;

class ARB_SYMBOL_VISIBLE network_label_dict;

class ARB_SYMBOL_VISIBLE network_selection;

class ARB_SYMBOL_VISIBLE network_value {
public:
    using custom_func_type = std::function<double(const network_connection_info&)>;

    network_value() { *this = network_value::scalar(0.0); }

    // Scalar value with conversion from double
    network_value(double value) { *this = network_value::scalar(value); }

    // Scalar value. Will always return the same value given at construction.
    static network_value scalar(double value);

    // A named value inside a network label dictionary
    static network_value named(std::string name);

    // Distamce netweem source and destination site
    static network_value distance(double scale = 1.0);

    // Uniform random value in (range[0], range[1]].
    static network_value uniform_distribution(unsigned seed, const std::array<double, 2>& range);

    // Radom value from a normal distribution with given mean and standard deviation.
    static network_value normal_distribution(unsigned seed, double mean, double std_deviation);

    // Radom value from a truncated normal distribution with given mean and standard deviation (of a
    // non-truncated normal distribution), where the value is always in (range[0], range[1]].
    // Note: Values are generated by reject-accept method from a normal
    // distribution. Low acceptance rate can leed to poor performance, for example with very small
    // ranges or a mean far outside the range.
    static network_value truncated_normal_distribution(unsigned seed,
        double mean,
        double std_deviation,
        const std::array<double, 2>& range);

    // Custom value using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result.
    static network_value custom(custom_func_type func);

    static network_value add(network_value left, network_value right);

    static network_value sub(network_value left, network_value right);

    static network_value mul(network_value left, network_value right);

    static network_value div(network_value left, network_value right);

    static network_value exp(network_value v);

    static network_value log(network_value v);

    static network_value min(network_value left, network_value right);

    static network_value max(network_value left, network_value right);

    // if contained in selection, the true_value is used and the false_value otherwise.
    static network_value if_else(network_selection cond,
        network_value true_value,
        network_value false_value);

    ARB_ARBOR_API friend std::ostream& operator<<(std::ostream& os, const network_value& v);

private:
    network_value(std::shared_ptr<network_value_impl> impl);

    friend std::shared_ptr<network_value_impl> thingify(network_value v,
        const network_label_dict& dict);

    std::shared_ptr<network_value_impl> impl_;
};

ARB_ARBOR_API inline network_value operator+(network_value a, network_value b) {
    return network_value::add(std::move(a), std::move(b));
}

ARB_ARBOR_API inline network_value operator-(network_value a, network_value b) {
    return network_value::sub(std::move(a), std::move(b));
}

ARB_ARBOR_API inline network_value operator*(network_value a, network_value b) {
    return network_value::mul(std::move(a), std::move(b));
}

ARB_ARBOR_API inline network_value operator/(network_value a, network_value b) {
    return network_value::div(std::move(a), std::move(b));
}

ARB_ARBOR_API inline network_value operator+(network_value a) { return a; }

ARB_ARBOR_API inline network_value operator-(network_value a) {
    return network_value::mul(-1.0, std::move(a));
}

class ARB_SYMBOL_VISIBLE network_selection {
public:
    using custom_func_type = std::function<bool(const network_connection_info&)>;

    network_selection() { *this = network_selection::none(); }

    // Select all
    static network_selection all();

    // Select none
    static network_selection none();

    // Named selection in the network label dictionary
    static network_selection named(std::string name);

    // Only select connections between different cells
    static network_selection inter_cell();

    // Select connections with the given source cell kind
    static network_selection source_cell_kind(cell_kind kind);

    // Select connections with the given destination cell kind
    static network_selection destination_cell_kind(cell_kind kind);

    // Select connections with the given source label
    static network_selection source_label(std::vector<cell_tag_type> labels);

    // Select connections with the given destination label
    static network_selection destination_label(std::vector<cell_tag_type> labels);

    // Select connections with source cells matching the indices in the list
    static network_selection source_cell(std::vector<cell_gid_type> gids);

    // Select connections with source cells matching the indices in the range
    static network_selection source_cell(gid_range range);

    // Select connections with destination cells matching the indices in the list
    static network_selection destination_cell(std::vector<cell_gid_type> gids);

    // Select connections with destination cells matching the indices in the range
    static network_selection destination_cell(gid_range range);

    // Select connections that form a chain, such that source cell "i" is connected to the destination cell "i+1"
    static network_selection chain(std::vector<cell_gid_type> gids);

    // Select connections that form a chain, such that source cell "i" is connected to the destination cell "i+1"
    static network_selection chain(gid_range range);

    // Select connections that form a reversed chain, such that source cell "i+1" is connected to the destination cell "i"
    static network_selection chain_reverse(gid_range range);

    // Select connections, that are selected by both "left" and "right"
    static network_selection intersect(network_selection left, network_selection right);

    // Select connections, that are selected by either or both "left" and "right"
    static network_selection join(network_selection left, network_selection right);

    // Select connections, that are selected by "left", unless selected by "right"
    static network_selection difference(network_selection left, network_selection right);

    // Select connections, that are selected by "left" or "right", but not both
    static network_selection symmetric_difference(network_selection left, network_selection right);

    // Invert the selection
    static network_selection complement(network_selection s);

    // Random selection using the bernoulli random distribution with probability "p" between 0.0
    // and 1.0
    static network_selection random(unsigned seed, network_value p);

    // Custom selection using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result. For gap junction selection,
    // "func" must be symmetric (func(a,b) = func(b,a)).
    static network_selection custom(custom_func_type func);

    // only select within given distance. This may enable more efficient sampling through an
    // internal spatial data structure.
    static network_selection distance_lt(double d);

    // only select if distance greater then given distance. This may enable more efficient sampling
    // through an internal spatial data structure.
    static network_selection distance_gt(double d);

    ARB_ARBOR_API friend std::ostream& operator<<(std::ostream& os, const network_selection& s);

private:
    network_selection(std::shared_ptr<network_selection_impl> impl);

    friend std::shared_ptr<network_selection_impl> thingify(network_selection s,
        const network_label_dict& dict);

    friend class network_value;

    std::shared_ptr<network_selection_impl> impl_;
};

class ARB_SYMBOL_VISIBLE network_label_dict {
public:
    using ns_map = std::unordered_map<std::string, network_selection>;
    using nv_map = std::unordered_map<std::string, network_value>;

    // Store a network selection under the given name
    network_label_dict& set(const std::string& name, network_selection s);

    // Store a network value under the given name
    network_label_dict& set(const std::string& name, network_value v);

    // Returns the stored network selection of the given name if it exists. None otherwise.
    std::optional<network_selection> selection(const std::string& name) const;

    // Returns the stored network value of the given name if it exists. None otherwise.
    std::optional<network_value> value(const std::string& name) const;

    // All stored network selections
    inline const ns_map& selections() const { return selections_; }

    // All stored network value
    inline const nv_map& values() const { return values_; }

private:
    ns_map selections_;
    nv_map values_;
};

// A complete network description required for processing
struct ARB_SYMBOL_VISIBLE network_description {
    network_selection selection;
    network_value weight;
    network_value delay;
    network_label_dict dict;
};

// Join two network selections
ARB_ARBOR_API network_selection join(network_selection left, network_selection right);

// Join three or more network selections
ARB_ARBOR_API network_selection join(network_selection left, network_selection right);
template <typename... Args>
network_selection join(network_selection l, network_selection r, Args... args) {
    return join(join(std::move(l), std::move(r)), std::move(args)...);
}

// Intersect two network selections
ARB_ARBOR_API network_selection join(network_selection left, network_selection right);
ARB_ARBOR_API network_selection intersect(network_selection left, network_selection right);

// Intersect three or more network selections
ARB_ARBOR_API network_selection join(network_selection left, network_selection right);
template <typename... Args>
network_selection intersect(network_selection l, network_selection r, Args... args) {
    return intersect(intersect(std::move(l), std::move(r)), std::move(args)...);
}

}  // namespace arb
