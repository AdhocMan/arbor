#pragma once

#include <arbor/common_types.hpp>
#include <arbor/export.hpp>
#include <arbor/recipe.hpp>

#include <optional>
#include <array>
#include <functional>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <utility>
#include <tuple>

namespace arb {

using network_location = std::array<double, 3>;

struct ARB_SYMBOL_VISIBLE network_cell_group {

    network_cell_group(cell_gid_type gid_begin,
        cell_gid_type gid_end,
        std::vector<cell_local_label_type> src_labels,
        std::vector<cell_local_label_type> dest_labels);

    cell_gid_type gid_begin, gid_end;
    std::vector<cell_local_label_type> src_labels, dest_labels;
};

struct ARB_SYMBOL_VISIBLE spatial_network_cell_group {
    spatial_network_cell_group(network_cell_group group, std::vector<network_location> locations);

    spatial_network_cell_group(cell_gid_type gid_begin,
        std::vector<cell_local_label_type> src_labels,
        std::vector<cell_local_label_type> dest_labels,
        std::vector<network_location> locations);

    inline const network_location& location(cell_gid_type gid) const {
        return locations.at(gid - group.gid_begin);
    }

    network_cell_group group;
    std::vector<network_location> locations;
};

struct ARB_SYMBOL_VISIBLE network_gj_group {

    network_gj_group(cell_gid_type gid_begin,
        cell_gid_type gid_end,
        std::vector<cell_local_label_type> labels);

    cell_gid_type gid_begin, gid_end;
    std::vector<cell_local_label_type> labels;
};

struct ARB_SYMBOL_VISIBLE spatial_network_gj_group {
    spatial_network_gj_group(network_gj_group group, std::vector<network_location> locations);

    spatial_network_gj_group(cell_gid_type gid_begin,
        std::vector<cell_local_label_type> labels,
        std::vector<network_location> locations);

    inline const network_location& location(cell_gid_type gid) const {
        return locations.at(gid - group.gid_begin);
    }

    network_gj_group group;
    std::vector<network_location> locations;
};

template <typename T,
    typename = std::enable_if_t<std::is_same_v<T, network_cell_group> ||
                                std::is_same_v<T, network_gj_group>>>
struct spatial {
    spatial(T group, std::vector<network_location> locations):
        group(std::move(group)),
        locations(std::move(locations)) {
        if (this->locations.size() != this->group.gid_end - this->group.gid_begin)
            throw std::runtime_error("spatial network group: The number of points is not "
                                     "equal to the network group size.");
    }

    inline const network_location& location(cell_gid_type gid) const {
        return locations.at(gid - group.gid_begin);
    }

    T group;
    std::vector<network_location> locations;
};

class ARB_SYMBOL_VISIBLE spatial_network_selection;

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
    bool operator()(const cell_global_label_type& src, const cell_global_label_type& dest) const;

private:
    friend cell_connection_network;
    friend gap_junction_network;
    friend spatial_network_selection;

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


class ARB_SYMBOL_VISIBLE spatial_network_selection {
public:
    spatial_network_selection(network_selection s);

    // Custom selection using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result. For gap junction selection,
    // "func" must be symmetric (func(a,b) = func(b,a)).
    static spatial_network_selection custom(std::function<bool(const cell_global_label_type&,
            const network_location&,
            const cell_global_label_type&,
            const network_location&,
            double)> func);

    static spatial_network_selection within_distance(double distance);

    // Returns true if a connection between src and dest is selected.
    bool operator()(const cell_global_label_type& src,
        const network_location& src_location,
        const cell_global_label_type& dest,
        const network_location& dest_location) const;

    spatial_network_selection operator&(spatial_network_selection right) const;

    spatial_network_selection operator|(spatial_network_selection right) const;

    spatial_network_selection operator^(spatial_network_selection right) const;

private:
    friend cell_connection_network;
    friend gap_junction_network;

    struct spatial_selection_impl {
        virtual bool select(cell_gid_type,
            const cell_local_label_type&,
            const network_location&,
            cell_gid_type,
            const cell_local_label_type&,
            const network_location&, double distance) const = 0;

        virtual std::optional<double> max_distance() const { return std::nullopt; }

        virtual ~spatial_selection_impl() = default;
    };

    struct selection_conversion_impl;
    struct and_impl;
    struct or_impl;
    struct xor_impl;
    struct custom_impl;
    struct within_distance_impl;

    spatial_network_selection(std::shared_ptr<spatial_selection_impl> impl);

    inline std::optional<double> max_distance() const { return impl_->max_distance(); }

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location, double distance) const {
        return impl_->select(
            src_gid, src_label, src_location, dest_gid, dest_label, dest_location, distance);
    }

    std::shared_ptr<spatial_selection_impl> impl_;
};

inline spatial_network_selection operator&(const network_selection& lhs,
    const spatial_network_selection& rhs) {
    return spatial_network_selection(lhs) & rhs;
}

inline spatial_network_selection operator&(const spatial_network_selection& lhs,
    const network_selection& rhs) {
    return lhs & spatial_network_selection(rhs);
}

inline spatial_network_selection operator|(const network_selection& lhs,
    const spatial_network_selection& rhs) {
    return spatial_network_selection(lhs) | rhs;
}

inline spatial_network_selection operator|(const spatial_network_selection& lhs,
    const network_selection& rhs) {
    return lhs | spatial_network_selection(rhs);
}

inline spatial_network_selection operator^(const network_selection& lhs,
    const spatial_network_selection& rhs) {
    return spatial_network_selection(lhs) ^ rhs;
}

inline spatial_network_selection operator^(const spatial_network_selection& lhs,
    const network_selection& rhs) {
    return lhs ^ spatial_network_selection(rhs);
}

class ARB_SYMBOL_VISIBLE spatial_network_value;

class ARB_SYMBOL_VISIBLE network_value {
public:
    // Uniform value
    network_value(double value);

    // Uniform random value in (range[0], range[1]]. Always returns the same value for repeated
    // calls with the same arguments and calls are symmetric v(a, b) = v(b, a).
    static network_value uniform_distribution(unsigned seed, const std::array<double, 2>& range);

    // Radom value from a normal distribution with given mean and standard deviation. Always returns
    // the same value for repeated calls with the same arguments and calls are symmetric v(a, b) =
    // v(b, a).
    static network_value normal_distribution(unsigned seed, double mean, double std_deviation);

    // Radom value from a truncated normal distribution with given mean and standard deviation (of a
    // non-truncated normal distribution), where the value is always in (range[0], range[1]]. Always
    // returns the same value for repeated calls with the same arguments and calls are symmetric
    // v(a, b) = v(b, a). Note: Values are generated by reject-accept method from a normal
    // distribution. Low acceptance rate can leed to poor performance, for example with very small
    // ranges or a mean far outside the range.
    static network_value truncated_normal_distribution(unsigned seed,
        double mean,
        double std_deviation,
        const std::array<double, 2>& range);

    // Custom value using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result. For gap junction values,
    // "func" must be symmetric (func(a,b) = func(b,a)).
    static network_value custom(
        std::function<double(const cell_global_label_type&, const cell_global_label_type&)> func);

    // Uniform value. Will always return the same value given at construction.
    static network_value uniform(double value);

    // Returns the value for a connection between src and dest.
    inline double operator()(const cell_global_label_type& src,
        const cell_global_label_type& dest) const {
        return impl_->get(src.gid, src.label, dest.gid, dest.label);
    }

private:
    friend cell_connection_network;
    friend gap_junction_network;
    friend spatial_network_value;

    struct value_impl {
        virtual double get(cell_gid_type,
            const cell_local_label_type&,
            cell_gid_type,
            const cell_local_label_type&) const = 0;

        virtual ~value_impl() = default;
    };

    struct uniform_distribution_impl;
    struct normal_distribution_impl;
    struct truncated_normal_distribution_impl;
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

class ARB_SYMBOL_VISIBLE spatial_network_value {
public:
    spatial_network_value(network_value value);

    // Custom value using the provided function "func". Repeated calls with the same arguments
    // to "func" must yield the same result. For gap junction values,
    // "func" must be symmetric (func(a,b) = func(b,a)).
    static spatial_network_value custom(std::function<double(const cell_global_label_type&,
            const network_location&,
            const cell_global_label_type&,
            const network_location&,
            double)> func);

    // Returns the value for a connection between src and dest.
    inline double operator()(const cell_global_label_type& src,
        const network_location& src_location,
        const cell_global_label_type& dest,
        const network_location& dest_location) const;

private:
    friend cell_connection_network;
    friend gap_junction_network;

    struct spatial_value_impl {
        virtual double get(cell_gid_type,
            const cell_local_label_type&,
            const network_location&,
            cell_gid_type,
            const cell_local_label_type&,
            const network_location&,
            double distance) const = 0;

        virtual ~spatial_value_impl() = default;
    };

    struct custom_impl;
    struct value_conversion_impl;

    spatial_network_value(std::shared_ptr<spatial_value_impl> impl);

    inline double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location, double distance) const {
        return impl_->get(
            src_gid, src_label, src_location, dest_gid, dest_label, dest_location, distance);
    }

    std::shared_ptr<spatial_value_impl> impl_;
};


// Generate cell connections on demand
class ARB_SYMBOL_VISIBLE cell_connection_network {
public:
    cell_connection_network();

    cell_connection_network(network_value weight,
        network_value delay,
        network_selection selection,
        std::vector<network_cell_group> pop);

    cell_connection_network(spatial_network_value weight,
        spatial_network_value delay,
        spatial_network_selection selection,
        std::vector<spatial_network_cell_group> pop);

    // Generate connections for the given global cell index
    inline std::vector<cell_connection> generate(cell_gid_type gid) const {
        return impl_->generate(gid);
    }

private:
    struct cell_connection_network_impl {
        virtual std::vector<cell_connection> generate(cell_gid_type gid) const = 0;

        virtual ~cell_connection_network_impl() = default;
    };

    struct empty_impl;
    struct spatial_impl;
    struct non_spatial_impl;

    std::shared_ptr<cell_connection_network_impl> impl_;
};


// Generate gap junctions on demand
class ARB_SYMBOL_VISIBLE gap_junction_network {
public:
    gap_junction_network();

    gap_junction_network(network_value weight,
        network_selection selection,
        std::vector<network_gj_group> pop);

    gap_junction_network(spatial_network_value weight,
        spatial_network_selection selection,
        std::vector<spatial_network_gj_group> pop);

    // Generate connections for the given global cell index
    inline std::vector<gap_junction_connection> generate(cell_gid_type gid) const {
        return impl_->generate(gid);
    }

private:
    struct gap_junction_network_impl {
        virtual std::vector<gap_junction_connection> generate(cell_gid_type gid) const = 0;

        virtual ~gap_junction_network_impl() = default;
    };

    struct empty_impl;
    struct spatial_impl;
    struct non_spatial_impl;

    std::shared_ptr<gap_junction_network_impl> impl_;
};

}  // namespace arb
