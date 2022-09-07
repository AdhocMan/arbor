#pragma once

#include <arbor/common_types.hpp>
#include <arbor/export.hpp>
#include <arbor/recipe.hpp>

#include <array>
#include <functional>
#include <memory>
#include <ostream>
#include <vector>

namespace arb {

using network_population = std::vector<cell_global_range_label_type>;

// Removes any duplicate global labels contained by merging overlapping intervals for matching local
// labels.
ARB_ARBOR_API network_population unique(const network_population& pop);

// Combine network populations. Does not remove duplicates.
ARB_ARBOR_API network_population join(const network_population& a, network_population b);

template <typename... ARGS>
inline network_population join(const network_population& a, network_population b, ARGS&&... args) {
    return join(join(a, b), std::forward<ARGS>(args)...);
}

class ARB_SYMBOL_VISIBLE cell_connection_network;

class ARB_SYMBOL_VISIBLE gap_junction_network;

class ARB_SYMBOL_VISIBLE network_selection {
public:
    // Random selection using the bernoulli random distribution with probability "p" between 0.0
    // and 1.0
    static network_selection bernoulli_random(unsigned seed, double p);

    // Custom selection using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result. For gap junction selection,
    // "func" must be symmetric (func(a,b) = func(b,a)).
    static network_selection custom(
        std::function<bool(const cell_global_label_type&, const cell_global_label_type&)> func);

    // Select all
    static network_selection all();

    // Select none
    static network_selection none();

    // Invert the selection
    static network_selection invert(network_selection s);

    // Only select connections between different cells
    static network_selection inter_cell();

    // Only select connections when the global labels are not equal. May select intra-cell
    // connections, if the local label is not equal.
    static network_selection not_equal();

    network_selection operator&(network_selection right) const;

    network_selection operator|(network_selection right) const;

    network_selection operator^(network_selection right) const;

    // Returns true if a connection between src and dest is selected.
    inline bool operator()(
        const cell_global_label_type& src,
        const cell_global_label_type& dest) const {
        return impl_->select(src.gid, src.label, dest.gid, dest.label);
    }

private:
    friend cell_connection_network;
    friend gap_junction_network;

    struct selection_impl {
        virtual bool select(cell_gid_type,
            const cell_local_label_type&,
            cell_gid_type,
            const cell_local_label_type&) const = 0;

        virtual ~selection_impl() = default;
    };

    struct bernoulli_random_impl;
    struct custom_impl;
    struct inter_cell_impl;
    struct not_equal_impl;
    struct all_impl;
    struct none_impl;
    struct and_impl;
    struct or_impl;
    struct xor_impl;
    struct invert_impl;

    network_selection(std::shared_ptr<selection_impl> impl);

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const {
        return impl_->select(src_gid, src_label, dest_gid, dest_label);
    }

    std::shared_ptr<selection_impl> impl_;
};

class ARB_SYMBOL_VISIBLE network_value {
public:
    // Uniform value
    network_value(double value);

    // Uniform random value in (range[0], range[1]]. Always returns the same value for repeated
    // calls with the same arguments and calls are symmetric v(a, b) = v(b, a).
    static network_value uniform_random(unsigned seed, std::array<double, 2> range);

    // Custom value using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result. For gap junction values,
    // "func" must be symmetric (func(a,b) = func(b,a)).
    static network_value custom(
        std::function<double(const cell_global_label_type&, const cell_global_label_type&)> func);

    // Uniform value. Will always return the same value given at construction.
    static network_value uniform(double value);

    // Returns the value for a connection between src and dest.
    inline double operator()(
        const cell_global_label_type& src,
        const cell_global_label_type& dest) const {
        return impl_->get(src.gid, src.label, dest.gid, dest.label);
    }

private:
    friend cell_connection_network;
    friend gap_junction_network;

    struct value_impl {
        virtual double get(cell_gid_type,
            const cell_local_label_type&,
            cell_gid_type,
            const cell_local_label_type&) const = 0;

        virtual ~value_impl() = default;
    };

    struct uniform_random_impl;
    struct custom_impl;
    struct uniform_impl;

    network_value(std::shared_ptr<value_impl> impl);

    inline double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const {
        return impl_->get(src_gid, src_label, dest_gid, dest_label);
    }

    std::shared_ptr<value_impl> impl_;
};

// Generate cell connections on demand
class ARB_SYMBOL_VISIBLE cell_connection_network {
public:
    cell_connection_network(network_value weight,
        network_value delay,
        network_selection selection,
        network_population src_pop,
        network_population dest_pop);

    inline const network_value& weight() const { return weight_; }

    inline const network_value& delay() const { return delay_; }

    inline const network_selection& selection() const { return selection_; }

    inline const network_population& source_population() const { return src_pop_; }

    inline const network_population& destination_population() const { return dest_pop_; }

    // Generate connections for the given global cell index
    std::vector<cell_connection> generate(cell_gid_type gid) const;

private:
    network_selection selection_;
    network_value weight_, delay_;
    network_population src_pop_;
    network_population dest_pop_;
};

// Generate gap junctions on demand
class ARB_SYMBOL_VISIBLE gap_junction_network {
public:
    gap_junction_network(network_value weight,
        network_selection selection,
        network_population src_pop,
        network_population dest_pop);

    inline const network_value& weight() const { return weight_; }

    inline const network_selection& selection() const { return selection_; }

    inline const network_population& source_population() const { return src_pop_; }

    inline const network_population& destination_population() const { return dest_pop_; }

    // Generate connections for the given global cell index
    std::vector<gap_junction_connection> generate(cell_gid_type gid) const;

private:
    network_selection selection_;
    network_value weight_;
    network_population src_pop_;
    network_population dest_pop_;
};


}  // namespace arb
