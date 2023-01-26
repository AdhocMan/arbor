#include <arbor/common_types.hpp>
#include <arbor/export.hpp>
#include <arbor/network.hpp>
#include <arbor/recipe.hpp>
#include <stdexcept>
#include <util/spatial_tree.hpp>

#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <cstdint>
#include <cmath>

#include "Random123/boxmuller.hpp"
#include "util/strprintf.hpp"

#include <Random123/threefry.h>
#include <Random123/uniform.hpp>

namespace arb {

namespace {

// Partial seed to use for network_value and network_selection generation.
// Different seed for each type to avoid unintentional correlation.
enum class network_seed : unsigned {
    selection_bernoulli = 2058443,
    value_uniform = 48202,
    value_normal = 8405,
    value_truncated_normal = 380237
};

// We only need minimal hash collisions and good spread over the hash range, because this will be
// used as input for random123, which then provides all desired hash properties.
// std::hash is implementation dependent, so we define our own for reproducibility.

std::uint64_t simple_string_hash(const std::string& s) {
    // use fnv1a hash algorithm
    constexpr std::uint64_t prime = 1099511628211ull;
    std::uint64_t h = 14695981039346656037ull;

    for (auto c: s) {
        h ^= c;
        h *= prime;
    }

    return h;
}

std::uint64_t combine_hash_pair(std::uint64_t a, std::uint64_t b) {
    // Use golden ration hashing
    constexpr std::uint64_t golden_ratio = 11400714819323198485llu;
    return a * golden_ratio + b;
}

std::uint64_t hash_global_tag(cell_gid_type gid, const std::string& tag) {
    return combine_hash_pair(gid, simple_string_hash(tag));
}

double uniform_rand_from_key_pair(std::array<unsigned, 2> seed,
    std::uint64_t key_a,
    std::uint64_t key_b) {
    using rand_type = r123::Threefry2x64;
    const rand_type::ctr_type seed_input = {{seed[0], seed[1]}};

    const rand_type::key_type key = {{std::min(key_a, key_b), std::max(key_a, key_b)}};
    rand_type gen;
    return r123::u01<double>(gen(seed_input, key)[0]);
}

double normal_rand_from_key_pair(std::array<unsigned, 2> seed,
    std::uint64_t key_a,
    std::uint64_t key_b) {
    using rand_type = r123::Threefry2x64;
    const rand_type::ctr_type seed_input = {{seed[0], seed[1]}};

    const rand_type::key_type key = {{std::min(key_a, key_b), std::max(key_a, key_b)}};
    rand_type gen;
    const auto rand_num = gen(seed_input, key);
    return r123::boxmuller(rand_num[0], rand_num[1]).x;
}

}  // namespace

network_cell_group::network_cell_group(cell_gid_type gid_begin,
    cell_gid_type gid_end,
    std::vector<cell_local_label_type> src_labels,
    std::vector<cell_local_label_type> dest_labels):
    gid_begin(gid_begin),
    gid_end(gid_end),
    src_labels(std::move(src_labels)),
    dest_labels(std::move(dest_labels)) {}

spatial_network_cell_group::spatial_network_cell_group(network_cell_group group,
    std::vector<network_location> locations):
    gid_begin(group.gid_begin),
    gid_end(group.gid_end),
    src_labels(std::move(group.src_labels)),
    dest_labels(std::move(group.dest_labels)),
    locations(std::move(locations)) {
    if (this->locations.size() != this->gid_end - this->gid_begin)
        throw std::runtime_error("spatial_network_cell_group: The number of points is not "
                                 "equal to the network cell group size.");
}

spatial_network_cell_group::spatial_network_cell_group(cell_gid_type gid_begin,
    std::vector<cell_local_label_type> src_labels,
    std::vector<cell_local_label_type> dest_labels,
    std::vector<network_location> locations):
    gid_begin(gid_begin),
    gid_end(gid_begin + locations.size()),
    src_labels(std::move(src_labels)),
    dest_labels(std::move(dest_labels)),
    locations(std::move(locations)) {}

network_gj_group::network_gj_group(cell_gid_type gid_begin,
    cell_gid_type gid_end,
    std::vector<cell_local_label_type> labels):
    gid_begin(gid_begin),
    gid_end(gid_end),
    labels(std::move(labels)) {}

spatial_network_gj_group::spatial_network_gj_group(network_gj_group group,
    std::vector<network_location> locations):
    gid_begin(group.gid_begin),
    gid_end(group.gid_end),
    labels(std::move(group.labels)),
    locations(std::move(locations)) {
    if (this->locations.size() != this->gid_end - this->gid_begin)
        throw std::runtime_error("spatial_network_cell_group: The number of points is not "
                                 "equal to the network cell group size.");
}

spatial_network_gj_group::spatial_network_gj_group(cell_gid_type gid_begin,
    std::vector<cell_local_label_type> labels,
    std::vector<network_location> locations):
    gid_begin(gid_begin),
    gid_end(gid_begin + locations.size()),
    labels(std::move(labels)),
    locations(std::move(locations)) {}

struct network_selection::bernoulli_random_impl: public selection_impl {
    double p = 0;
    unsigned seed = 0;

    bernoulli_random_impl(double probability, unsigned rand_seed):
        p(probability),
        seed(rand_seed) {}

    bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return uniform_rand_from_key_pair({unsigned(network_seed::selection_bernoulli), seed},
                   hash_global_tag(src_gid, src_label.tag),
                   hash_global_tag(dest_gid, dest_label.tag)) < p;
    }
};

struct network_selection::inter_cell_impl: public selection_impl {
    bool select(cell_gid_type src_gid,
        const cell_local_label_type&,
        cell_gid_type dest_gid,
        const cell_local_label_type&) const override {
        return src_gid != dest_gid;
    }
};

struct network_selection::not_equal_impl: public selection_impl {
    bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return (src_gid != dest_gid) || (src_label != dest_label);
    }
};

struct network_selection::all_impl: public selection_impl {
    bool select(cell_gid_type,
        const cell_local_label_type&,
        cell_gid_type,
        const cell_local_label_type&) const override {
        return true;
    }
};

struct network_selection::none_impl: public selection_impl {
    bool select(cell_gid_type,
        const cell_local_label_type&,
        cell_gid_type,
        const cell_local_label_type&) const override {
        return false;
    }
};

struct network_selection::custom_impl: public selection_impl {
    std::function<bool(const cell_global_label_type&, const cell_global_label_type&)> func;

    custom_impl(
        std::function<bool(const cell_global_label_type&, const cell_global_label_type&)> f):
        func(std::move(f)) {}

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return func({src_gid, src_label}, {dest_gid, dest_label});
    }
};

struct network_selection::and_impl: public selection_impl {
    network_selection left, right;

    and_impl(network_selection l, network_selection r): left(std::move(l)), right(std::move(r)) {}


    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return left.select(src_gid, src_label, dest_gid, dest_label) &&
               right.select(src_gid, src_label, dest_gid, dest_label);
    }
};

struct network_selection::or_impl: public selection_impl {
    network_selection left, right;

    or_impl(network_selection l, network_selection r): left(std::move(l)), right(std::move(r)) {}

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return left.select(src_gid,
                   src_label,
                   dest_gid,
                   dest_label) ||
               right.select(
                   src_gid, src_label, dest_gid, dest_label);
    }
};

struct network_selection::xor_impl: public selection_impl {
    network_selection left, right;

    xor_impl(network_selection l, network_selection r): left(std::move(l)), right(std::move(r)) {}

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return left.select(src_gid, src_label, dest_gid, dest_label) ^
               right.select(src_gid, src_label, dest_gid, dest_label);
    }
};

struct network_selection::invert_impl: public selection_impl {
    network_selection selection;

    invert_impl(network_selection s): selection(std::move(s)) {}

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return !(selection.select(src_gid, src_label, dest_gid, dest_label));
    }
};


network_selection::network_selection(std::shared_ptr<selection_impl> impl):
    impl_(std::move(impl)) {}

network_selection network_selection::bernoulli_random(unsigned seed, double p) {
    return {std::shared_ptr<selection_impl>(new bernoulli_random_impl(p, seed))};
}

network_selection network_selection::custom(
    std::function<bool(const cell_global_label_type&, const cell_global_label_type&)> func) {
    return {std::shared_ptr<selection_impl>(new network_selection::custom_impl(std::move(func)))};
}

network_selection network_selection::inter_cell() {
    return {std::shared_ptr<selection_impl>(new network_selection::inter_cell_impl())};
}

network_selection network_selection::not_equal() {
    return {std::shared_ptr<selection_impl>(new network_selection::not_equal_impl())};
}

network_selection network_selection::all() {
    return {std::shared_ptr<selection_impl>(new network_selection::all_impl())};
}

network_selection network_selection::none() {
    return {std::shared_ptr<selection_impl>(new network_selection::none_impl())};
}

network_selection network_selection::invert(network_selection s) {
    return {std::shared_ptr<selection_impl>(new network_selection::invert_impl(std::move(s)))};
}

network_selection network_selection::operator&(network_selection right) const {
    return {
        std::shared_ptr<selection_impl>(new network_selection::and_impl(*this, std::move(right)))};
}

network_selection network_selection::operator|(network_selection right) const {
    return {
        std::shared_ptr<selection_impl>(new network_selection::or_impl(*this, std::move(right)))};
}

network_selection network_selection::operator^(network_selection right) const {
    return {
        std::shared_ptr<selection_impl>(new network_selection::xor_impl(*this, std::move(right)))};
}

bool network_selection::operator()(const cell_global_label_type& src,
    const cell_global_label_type& dest) const {
    return impl_->select(src.gid, src.label, dest.gid, dest.label);
}


struct spatial_network_selection::selection_conversion_impl: public spatial_selection_impl {
    network_selection selection;

    selection_conversion_impl(network_selection s): selection(std::move(s)) {}

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location&,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location&,
        double) const override {
        return selection.select(src_gid, src_label, dest_gid, dest_label);
    }
};

struct spatial_network_selection::and_impl: public spatial_selection_impl {
    spatial_network_selection left, right;

    and_impl(spatial_network_selection l, spatial_network_selection r): left(std::move(l)), right(std::move(r)) {}

    std::optional<double> max_distance() const override {
        const auto d_left = left.max_distance();
        const auto d_right = right.max_distance();

        if (d_left && d_right) return std::min(d_left.value(), d_right.value());
        if (d_left) return d_left.value();
        if (d_right) return d_right.value();

        return std::nullopt;
    }

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location,
        double distance) const override {
        return left.select(src_gid,
                   src_label,
                   src_location,
                   dest_gid,
                   dest_label,
                   dest_location,
                   distance) &&
               right.select(
                   src_gid, src_label, src_location, dest_gid, dest_label, dest_location, distance);
    }
};

struct spatial_network_selection::or_impl: public spatial_selection_impl {
    spatial_network_selection left, right;

    or_impl(spatial_network_selection l, spatial_network_selection r): left(std::move(l)), right(std::move(r)) {}

    std::optional<double> max_distance() const override {
        const auto d_left = left.max_distance();
        const auto d_right = right.max_distance();

        if (d_left && d_right) return std::max(d_left.value(), d_right.value());

        return std::nullopt;
    }

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location,
        double location) const override {
        return left.select(src_gid,
                   src_label,
                   src_location,
                   dest_gid,
                   dest_label,
                   dest_location,
                   location) ||
               right.select(
                   src_gid, src_label, src_location, dest_gid, dest_label, dest_location, location);
    }
};

struct spatial_network_selection::xor_impl: public spatial_selection_impl {
    spatial_network_selection left, right;

    xor_impl(spatial_network_selection l, spatial_network_selection r): left(std::move(l)), right(std::move(r)) {}

    std::optional<double> max_distance() const override {
        const auto d_left = left.max_distance();
        const auto d_right = right.max_distance();

        if (d_left && d_right) return std::max(d_left.value(), d_right.value());

        return std::nullopt;
    }

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location,
        double distance) const override {
        return left.select(src_gid,
                   src_label,
                   src_location,
                   dest_gid,
                   dest_label,
                   dest_location,
                   distance) ^
               right.select(
                   src_gid, src_label, src_location, dest_gid, dest_label, dest_location, distance);
    }
};


struct spatial_network_selection::custom_impl: public spatial_selection_impl {
    std::function<bool(const cell_global_label_type&,
        const network_location&,
        const cell_global_label_type&,
        const network_location&,
        double)>
        func;

    custom_impl(std::function<bool(const cell_global_label_type&,
            const network_location&,
            const cell_global_label_type&,
            const network_location&,
            double)> f):
        func(std::move(f)) {}

    inline bool select(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location,
        double distance) const override {
        return func({src_gid, src_label}, src_location, {dest_gid, dest_label}, dest_location, distance);
    }
};


struct spatial_network_selection::within_distance_impl: public spatial_selection_impl {
    double d;

    within_distance_impl(double d): d(d) {}

    std::optional<double> max_distance() const override { return d; }

    inline bool select(cell_gid_type,
        const cell_local_label_type&,
        const network_location&,
        cell_gid_type,
        const cell_local_label_type&,
        const network_location&,
        double distance) const override {
        return distance <= d;
    }
};

spatial_network_selection::spatial_network_selection(network_selection s):
    impl_(new selection_conversion_impl(std::move(s))) {}

spatial_network_selection::spatial_network_selection(std::shared_ptr<spatial_selection_impl> impl):
    impl_(std::move(impl)) {}

bool spatial_network_selection::operator()(const cell_global_label_type& src,
    const network_location& src_location,
    const cell_global_label_type& dest,
    const network_location& dest_location) const {
    const network_location diff = {src_location[0] - dest_location[0],
        src_location[1] - dest_location[1],
        src_location[2] - dest_location[2]};
    const double distance = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
    return impl_->select(src.gid, src.label, src_location, dest.gid, dest.label, dest_location, distance);
}

spatial_network_selection spatial_network_selection::operator&(
    spatial_network_selection right) const {
    return {std::shared_ptr<spatial_selection_impl>(
        new spatial_network_selection::and_impl(*this, std::move(right)))};
}

spatial_network_selection spatial_network_selection::operator|(
    spatial_network_selection right) const {
    return {std::shared_ptr<spatial_selection_impl>(
        new spatial_network_selection::or_impl(*this, std::move(right)))};
}

spatial_network_selection spatial_network_selection::operator^(
    spatial_network_selection right) const {
    return {std::shared_ptr<spatial_selection_impl>(
        new spatial_network_selection::xor_impl(*this, std::move(right)))};
}

spatial_network_selection spatial_network_selection::custom(
    std::function<bool(const cell_global_label_type&,
        const network_location&,
        const cell_global_label_type&,
        const network_location&,
        double)> func) {
    return {std::shared_ptr<spatial_selection_impl>(
        new spatial_network_selection::custom_impl(std::move(func)))};
}

spatial_network_selection spatial_network_selection::within_distance(double distance) {
    return {std::shared_ptr<spatial_selection_impl>(
        new spatial_network_selection::within_distance_impl(distance))};
}

struct network_value::uniform_distribution_impl: public value_impl {
    unsigned seed = 0;
    std::array<double, 2> range;

    uniform_distribution_impl(unsigned rand_seed, const std::array<double, 2>& r):
        seed(rand_seed),
        range(r) {
        if (range[0] >= range[1])
            throw std::invalid_argument("Uniform distribution: invalid range");
    }

    double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        if (range[0] > range[1]) return range[1];

        // random number between 0 and 1
        const auto rand_num =
            uniform_rand_from_key_pair({unsigned(network_seed::value_uniform), seed},
                hash_global_tag(src_gid, src_label.tag),
                hash_global_tag(dest_gid, dest_label.tag));

        return (range[1] - range[0]) * rand_num + range[0];
    }
};


struct network_value::normal_distribution_impl: public value_impl {
    unsigned seed = 0;
    double mean = 0.0;
    double std_deviation = 1.0;

    normal_distribution_impl(unsigned rand_seed, double mean_, double std_deviation_):
        seed(rand_seed),
        mean(mean_),
        std_deviation(std_deviation_) {}

    double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return mean + std_deviation *
                          normal_rand_from_key_pair({unsigned(network_seed::value_normal), seed},
                              hash_global_tag(src_gid, src_label.tag),
                              hash_global_tag(dest_gid, dest_label.tag));
    }
};

struct network_value::truncated_normal_distribution_impl: public value_impl {
    unsigned seed = 0;
    double mean = 0.0;
    double std_deviation = 1.0;
    std::array<double, 2> range;

    truncated_normal_distribution_impl(unsigned rand_seed,
        double mean_,
        double std_deviation_,
        const std::array<double, 2>& range_):
        seed(rand_seed),
        mean(mean_),
        std_deviation(std_deviation_),
        range(range_) {
        if (range[0] >= range[1])
            throw std::invalid_argument("Truncated normal distribution: invalid range");
    }

    double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {

        const auto src_hash = hash_global_tag(src_gid, src_label.tag);
        auto dest_hash = hash_global_tag(dest_gid, dest_label.tag);

        double value = 0.0;

        do {
            value =
                mean + std_deviation * normal_rand_from_key_pair(
                                           {unsigned(network_seed::value_truncated_normal), seed},
                                           src_hash,
                                           dest_hash);
            ++dest_hash;
        } while (!(value > range[0] && value <= range[1]));

        return value;
    }
};

struct network_value::custom_impl: public value_impl {
    std::function<double(const cell_global_label_type&, const cell_global_label_type&)> func;

    custom_impl(
        std::function<double(const cell_global_label_type&, const cell_global_label_type&)> f):
        func(std::move(f)) {}

    inline double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        return func({src_gid, src_label}, {dest_gid, dest_label});
    }
};

struct network_value::uniform_impl: public value_impl {
    double value;

    uniform_impl(double v): value(v) {}

    inline double get(cell_gid_type,
        const cell_local_label_type&,
        cell_gid_type,
        const cell_local_label_type&) const override {
        return value;
    }
};

network_value::network_value(double value): network_value(network_value::uniform(value)) {}

network_value::network_value(std::shared_ptr<value_impl> impl): impl_(std::move(impl)) {}

network_value network_value::uniform_distribution(unsigned seed,
    const std::array<double, 2>& range) {
    return {std::shared_ptr<value_impl>(new uniform_distribution_impl(seed, range))};
}

network_value network_value::normal_distribution(unsigned seed, double mean, double std_deviation) {
    return {std::shared_ptr<value_impl>(new normal_distribution_impl(seed, mean, std_deviation))};
}

network_value network_value::truncated_normal_distribution(unsigned seed,
    double mean,
    double std_deviation,
    const std::array<double, 2>& range) {
    return {std::shared_ptr<value_impl>(
        new truncated_normal_distribution_impl(seed, mean, std_deviation, range))};
}

network_value network_value::custom(
    std::function<double(const cell_global_label_type&, const cell_global_label_type&)> func) {
    return {std::shared_ptr<value_impl>(new custom_impl(std::move(func)))};
}

network_value network_value::uniform(double value) {
    return {std::shared_ptr<value_impl>(new uniform_impl(value))};
}


struct spatial_network_value::value_conversion_impl: public spatial_value_impl {
    network_value value;

    value_conversion_impl(network_value v): value(std::move(v)) {}

    inline double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location&,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location&,
        double) const override {
        return value.get(src_gid, src_label, dest_gid, dest_label);
    }
};

struct spatial_network_value::custom_impl: public spatial_value_impl {
    std::function<double(const cell_global_label_type&,
        const network_location&,
        const cell_global_label_type&,
        const network_location&,
        double)>
        func;

    custom_impl(std::function<double(const cell_global_label_type&,
            const network_location&,
            const cell_global_label_type&,
            const network_location&,
            double)> f):
        func(std::move(f)) {}

    inline double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        const network_location& src_location,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label,
        const network_location& dest_location,
        double distance) const override {
        return func({src_gid, src_label}, src_location, {dest_gid, dest_label}, dest_location, distance);
    }
};

spatial_network_value::spatial_network_value(network_value value):
    impl_(new value_conversion_impl(std::move(value))) {}

spatial_network_value::spatial_network_value(std::shared_ptr<spatial_value_impl> impl):
    impl_(std::move(impl)) {}

double spatial_network_value::operator()(const cell_global_label_type& src,
    const network_location& src_location,
    const cell_global_label_type& dest,
    const network_location& dest_location) const {
    const network_location diff = {src_location[0] - dest_location[0],
        src_location[1] - dest_location[1],
        src_location[2] - dest_location[2]};
    const double distance = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
    return impl_->get(src.gid, src.label, src_location, dest.gid, dest.label, dest_location, distance);
}

spatial_network_value spatial_network_value::custom(
    std::function<double(const cell_global_label_type&,
        const network_location&,
        const cell_global_label_type&,
        const network_location&,
        double)> func) {
    return {std::shared_ptr<spatial_value_impl>(
        new spatial_network_value::custom_impl(std::move(func)))};
}

struct cell_connection_network::empty_impl: public cell_connection_network_impl {
    std::vector<cell_connection> generate(cell_gid_type gid) const override { return {}; }
};

struct cell_connection_network::non_spatial_impl: public cell_connection_network_impl {
    network_selection selection_;
    network_value weight_, delay_;
    std::vector<network_cell_group> sorted_pop_;

    static bool group_comparison(const network_cell_group& a, const network_cell_group& b) {
                    return a.gid_begin < b.gid_begin;
    }

    non_spatial_impl(network_value weight,
        network_value delay,
        network_selection selection,
        std::vector<network_cell_group> pop):
        selection_(std::move(selection)),
        weight_(std::move(weight)),
        delay_(std::move(delay)),
        sorted_pop_(std::move(pop)) {

        if (!sorted_pop_.empty()) {
            // sort groups by first gid index
            std::sort(sorted_pop_.begin(), sorted_pop_.end(), group_comparison);

            // Make sure there is no overlap between groups
            for (std::size_t i = 0; i < sorted_pop_.size() - 1; ++i) {
                if(sorted_pop_[i].gid_end > sorted_pop_[i + 1].gid_begin)
                    throw std::runtime_error("Network cell groups must not overlap.");
            }
        }
    }

    // Generate connections for the given global cell index
    std::vector<cell_connection> generate(cell_gid_type gid) const override {
        std::vector<cell_connection> connections;

        // find potential group
        const auto group_it = std::lower_bound(sorted_pop_.begin(),
            sorted_pop_.end(),
            gid,
            [](const network_cell_group& g, cell_gid_type id) { return g.gid_begin < id && g.gid_end <= id; });

        // return empty connections if gid is not in group
        if (group_it == sorted_pop_.end() || group_it->gid_end <= gid) return connections;

        const auto& dest_labels = group_it->dest_labels;

        for (const auto& dest_label: dest_labels) {
            for (const auto& src_group: sorted_pop_) {
                for (const auto& src_label: src_group.src_labels) {
                    for (auto src_gid = src_group.gid_begin; src_gid < src_group.gid_end;
                         ++src_gid) {
                        if (selection_.select(src_gid, src_label, gid, dest_label))
                            connections.emplace_back(cell_global_label_type{src_gid, src_label},
                                dest_label,
                                weight_.get(src_gid, src_label, gid, dest_label),
                                delay_.get(src_gid, src_label, gid, dest_label));
                    }
                }
            }
        }

        // remove any duplicates
        std::sort(connections.begin(), connections.end());
        connections.erase(std::unique(connections.begin(), connections.end()), connections.end());

        return connections;
    }
};


struct cell_connection_network::spatial_impl: public cell_connection_network_impl {
    spatial_network_selection selection_;
    spatial_network_value weight_, delay_;
    std::vector<spatial_network_cell_group> sorted_pop_;
    std::optional<spatial_tree<std::pair<std::size_t, cell_gid_type>, 3>> src_tree_;

    static bool group_comparison(const spatial_network_cell_group& a, const spatial_network_cell_group& b) {
        return a.gid_begin < b.gid_begin;
    }

    spatial_impl(spatial_network_value weight,
        spatial_network_value delay,
        spatial_network_selection selection,
        std::vector<spatial_network_cell_group> pop):
        selection_(std::move(selection)),
        weight_(std::move(weight)),
        delay_(std::move(delay)),
        sorted_pop_(std::move(pop)) {

        if (!sorted_pop_.empty()) {
            // sort groups by first gid index
            std::sort(sorted_pop_.begin(), sorted_pop_.end(), group_comparison);

            // Make sure there is no overlap between groups
            for (std::size_t i = 0; i < sorted_pop_.size() - 1; ++i) {
                if(sorted_pop_[i].gid_end > sorted_pop_[i + 1].gid_begin)
                    throw std::runtime_error("Network cell groups must not overlap.");
            }

            if (selection_.max_distance().has_value()) {

                std::vector<std::pair<network_location, std::pair<std::size_t, cell_gid_type>>>
                    tree_data;

                std::size_t group_idx = 0;
                for(const auto& g : sorted_pop_) {
                    // only add if source labels exist
                    if (!g.src_labels.empty()) {
                        for (auto gid = g.gid_begin; gid < g.gid_end; ++gid) {
                            tree_data.emplace_back(g.locations[gid - g.gid_begin],
                                std::make_pair(group_idx, gid));
                        }
                    }
                    ++group_idx;
                }

                // Create tree with maximum depth of 8 and a leaf size target of 2^8
                src_tree_ = spatial_tree<std::pair<std::size_t, cell_gid_type>, 3>(
                    8, std::max<std::size_t>(10, tree_data.size() / 256), std::move(tree_data));
            }
        }

    }

    // Generate connections for the given global cell index
    std::vector<cell_connection> generate(cell_gid_type gid) const override {
        std::vector<cell_connection> connections;

        // find potential group
        const auto group_it = std::lower_bound(sorted_pop_.begin(),
            sorted_pop_.end(),
            gid,
            [](const spatial_network_cell_group& g, cell_gid_type id) {
                return g.gid_begin < id && g.gid_end <= id;
            });

        // return empty connections if gid is not in group
        if (group_it == sorted_pop_.end() || group_it->gid_end <= gid) return connections;

        const auto& dest_labels = group_it->dest_labels;
        const auto& dest_loc = group_it->location(gid);

        auto add_connections = [&](cell_gid_type src_gid,
                                   const network_location& src_loc,
                                   const std::vector<cell_local_label_type>& src_labels,
                                   cell_gid_type dest_gid,
                                   const network_location& dest_loc,
                                   const std::vector<cell_local_label_type>& dest_label) {
            if(dest_labels.empty() || src_labels.empty()) return;

            const network_location diff = {
                src_loc[0] - dest_loc[0], src_loc[1] - dest_loc[1], src_loc[2] - dest_loc[2]};
            const double distance =
                std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

            for (const auto& dest_label: dest_labels) {
                for (const auto& src_label: src_labels) {
                    if (selection_.select(
                            src_gid, src_label, src_loc, gid, dest_label, dest_loc, distance))
                        connections.emplace_back(cell_global_label_type{src_gid, src_label},
                            dest_label,
                            weight_.get(
                                src_gid, src_label, src_loc, gid, dest_label, dest_loc, distance),
                            delay_.get(
                                src_gid, src_label, src_loc, gid, dest_label, dest_loc, distance));
                }
            }
        };

        if (src_tree_.has_value()) {
            const auto d = selection_.max_distance().value();

            src_tree_.value().bounding_box_for_each(
                network_location{dest_loc[0] - d, dest_loc[1] - d, dest_loc[2] - d},
                network_location{dest_loc[0] + d, dest_loc[1] + d, dest_loc[2] + d},
                [&](const network_location& src_loc,
                    const std::pair<std::size_t, cell_gid_type>& src_indices) {
                    const auto src_gid = src_indices.second;

                    add_connections(src_gid,
                        src_loc,
                        sorted_pop_[src_indices.first].src_labels,
                        gid,
                        dest_loc,
                        dest_labels);
                });
        }
        else {
            for (const auto& src_group: sorted_pop_) {
                for (auto src_gid = src_group.gid_begin; src_gid < src_group.gid_end;
                     ++src_gid) {

                    const auto& src_loc = src_group.location(src_gid);

                    add_connections(src_gid,
                        src_loc,
                        src_group.src_labels,
                        gid,
                        dest_loc,
                        dest_labels);
                }
            }
        }

        // remove any duplicates
        std::sort(connections.begin(), connections.end());
        connections.erase(std::unique(connections.begin(), connections.end()), connections.end());

        return connections;
    }
};


cell_connection_network::cell_connection_network(): impl_(new cell_connection_network::empty_impl()) {}

cell_connection_network::cell_connection_network(network_value weight,
    network_value delay,
    network_selection selection,
    std::vector<network_cell_group> pop):
    impl_(new cell_connection_network::non_spatial_impl(std::move(weight),
        std::move(delay),
        std::move(selection),
        std::move(pop))) {}

cell_connection_network::cell_connection_network(spatial_network_value weight,
    spatial_network_value delay,
    spatial_network_selection selection,
    std::vector<spatial_network_cell_group> pop):
    impl_(new cell_connection_network::spatial_impl(std::move(weight),
        std::move(delay),
        std::move(selection),
        std::move(pop))) {}

struct gap_junction_network::empty_impl: public gap_junction_network_impl {
    std::vector<gap_junction_connection> generate(cell_gid_type gid) const override { return {}; }
};

struct gap_junction_network::non_spatial_impl: public gap_junction_network_impl {
    network_selection selection_;
    network_value weight_;
    std::vector<network_gj_group> sorted_pop_;

    static bool group_comparison(const network_gj_group& a, const network_gj_group& b) {
                    return a.gid_begin < b.gid_begin;
    }

    non_spatial_impl(network_value weight,
        network_selection selection,
        std::vector<network_gj_group> pop):
        selection_(std::move(selection)),
        weight_(std::move(weight)),
        sorted_pop_(std::move(pop)) {

        if (!sorted_pop_.empty()) {
            // sort groups by first gid index
            std::sort(sorted_pop_.begin(), sorted_pop_.end(), group_comparison);

            // Make sure there is no overlap between groups
            for (std::size_t i = 0; i < sorted_pop_.size() - 1; ++i) {
                if(sorted_pop_[i].gid_end > sorted_pop_[i + 1].gid_begin)
                    throw std::runtime_error("Network cell groups must not overlap.");
            }
        }
    }

    // Generate connections for the given global cell index
    std::vector<gap_junction_connection> generate(cell_gid_type gid) const override {
        std::vector<gap_junction_connection> connections;

        // find potential group
        const auto group_it = std::lower_bound(sorted_pop_.begin(),
            sorted_pop_.end(),
            gid,
            [](const network_gj_group& g, cell_gid_type id) {
                return g.gid_begin < id && g.gid_end <= id;
            });

        // return empty connections if gid is not in group
        if (group_it == sorted_pop_.end() || group_it->gid_end <= gid) return connections;

        const auto& dest_labels = group_it->labels;

        for (const auto& dest_label: dest_labels) {
            for (const auto& src_group: sorted_pop_) {
                for (const auto& src_label: src_group.labels) {
                    for (auto src_gid = src_group.gid_begin; src_gid < src_group.gid_end;
                         ++src_gid) {
                        if (selection_.select(src_gid, src_label, gid, dest_label))
                            connections.emplace_back(cell_global_label_type{src_gid, src_label},
                                dest_label,
                                weight_.get(src_gid, src_label, gid, dest_label));
                    }
                }
            }
        }

        // remove any duplicates
        std::sort(connections.begin(), connections.end());
        connections.erase(std::unique(connections.begin(), connections.end()), connections.end());

        return connections;
    }
};


struct gap_junction_network::spatial_impl: public gap_junction_network_impl {
    spatial_network_selection selection_;
    spatial_network_value weight_;
    std::vector<spatial_network_gj_group> sorted_pop_;
    std::optional<spatial_tree<std::pair<std::size_t, cell_gid_type>, 3>> tree_;

    static bool group_comparison(const spatial_network_gj_group& a, const spatial_network_gj_group& b) {
        return a.gid_begin < b.gid_begin;
    }

    spatial_impl(spatial_network_value weight,
        spatial_network_selection selection,
        std::vector<spatial_network_gj_group> pop):
        selection_(std::move(selection)),
        weight_(std::move(weight)),
        sorted_pop_(std::move(pop)) {

        if (!sorted_pop_.empty()) {
            // sort groups by first gid index
            std::sort(sorted_pop_.begin(), sorted_pop_.end(), group_comparison);

            // Make sure there is no overlap between groups
            for (std::size_t i = 0; i < sorted_pop_.size() - 1; ++i) {
                if(sorted_pop_[i].gid_end > sorted_pop_[i + 1].gid_begin)
                    throw std::runtime_error("Network gj groups must not overlap.");
            }

            if (selection_.max_distance().has_value()) {

                std::vector<std::pair<network_location, std::pair<std::size_t, cell_gid_type>>>
                    tree_data;

                std::size_t group_idx = 0;
                for(const auto& g : sorted_pop_) {
                    // only add if source labels exist
                    if (!g.labels.empty()) {
                        for (auto gid = g.gid_begin; gid < g.gid_end; ++gid) {
                            tree_data.emplace_back(g.locations[gid - g.gid_begin],
                                std::make_pair(group_idx, gid));
                        }
                    }
                    ++group_idx;
                }

                // Create tree with maximum depth of 8 and a leaf size target of 2^8
                tree_ = spatial_tree<std::pair<std::size_t, cell_gid_type>, 3>(
                    8, std::max<std::size_t>(10, tree_data.size() / 256), std::move(tree_data));
            }
        }
    }

    // Generate connections for the given global cell index
    std::vector<gap_junction_connection> generate(cell_gid_type gid) const override {
        std::vector<gap_junction_connection> connections;

        // find potential group
        const auto group_it = std::lower_bound(sorted_pop_.begin(),
            sorted_pop_.end(),
            gid,
            [](const spatial_network_gj_group& g, const cell_gid_type& id) {
                return g.gid_begin < id && g.gid_end <= id;
            });

        // return empty connections if gid is not in group
        if (group_it == sorted_pop_.end() || group_it->gid_end <= gid) return connections;

        const auto& dest_labels = group_it->labels;
        const auto& dest_loc = group_it->location(gid);

        auto add_connections = [&](cell_gid_type src_gid,
                                   const network_location& src_loc,
                                   const std::vector<cell_local_label_type>& src_labels,
                                   cell_gid_type dest_gid,
                                   const network_location& dest_loc,
                                   const std::vector<cell_local_label_type>& dest_label) {
            if(dest_labels.empty() || src_labels.empty()) return;

            const network_location diff = {
                src_loc[0] - dest_loc[0], src_loc[1] - dest_loc[1], src_loc[2] - dest_loc[2]};
            const double distance =
                std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

            if (selection_.max_distance() && selection_.max_distance().value() < distance) return;

            for (const auto& dest_label: dest_labels) {
                for (const auto& src_label: src_labels) {
                    if (selection_.select(
                            src_gid, src_label, src_loc, gid, dest_label, dest_loc, distance))
                        connections.emplace_back(cell_global_label_type{src_gid, src_label},
                            dest_label,
                            weight_.get(
                                src_gid, src_label, src_loc, gid, dest_label, dest_loc, distance));
                }
            }
        };

        if (tree_.has_value()) {
            const auto d = selection_.max_distance().value();

            tree_.value().bounding_box_for_each(
                network_location{dest_loc[0] - d, dest_loc[1] - d, dest_loc[2] - d},
                network_location{dest_loc[0] + d, dest_loc[1] + d, dest_loc[2] + d},
                [&](const network_location& src_loc,
                    const std::pair<std::size_t, cell_gid_type>& src_indices) {
                    const auto src_gid = src_indices.second;

                    add_connections(src_gid,
                        src_loc,
                        sorted_pop_[src_indices.first].labels,
                        gid,
                        dest_loc,
                        dest_labels);
                });
        }
        else {
            for (const auto& src_group: sorted_pop_) {
                for (auto src_gid = src_group.gid_begin; src_gid < src_group.gid_end;
                     ++src_gid) {

                    const auto& src_loc = src_group.location(src_gid);

                    add_connections(
                        src_gid, src_loc, src_group.labels, gid, dest_loc, dest_labels);
                }
            }
        }

        // remove any duplicates
        std::sort(connections.begin(), connections.end());
        connections.erase(std::unique(connections.begin(), connections.end()), connections.end());

        return connections;
    }
};


gap_junction_network::gap_junction_network(): impl_(new gap_junction_network::empty_impl()) {}

gap_junction_network::gap_junction_network(network_value weight,
    network_selection selection,
    std::vector<network_gj_group> pop):
    impl_(new gap_junction_network::non_spatial_impl(std::move(weight),
        std::move(selection),
        std::move(pop))) {}

gap_junction_network::gap_junction_network(spatial_network_value weight,
    spatial_network_selection selection,
    std::vector<spatial_network_gj_group> pop):
    impl_(new gap_junction_network::spatial_impl(std::move(weight),
        std::move(selection),
        std::move(pop))) {}

}  // namespace arb
