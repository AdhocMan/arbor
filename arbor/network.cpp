#include <arbor/common_types.hpp>
#include <arbor/export.hpp>
#include <arbor/network.hpp>
#include <arbor/recipe.hpp>

#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <cstdint>

#include "util/strprintf.hpp"

#include <Random123/threefry.h>
#include <Random123/uniform.hpp>

namespace arb {

namespace {

// Partial seed to use for network_value and network_selection generation.
// Different seed for each type to avoid unintentional correlation.
enum class network_seed : unsigned { value = 48202, selection = 2058443 };

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
}  // namespace

ARB_ARBOR_API network_population unique(const network_population& pop) {
    std::map<cell_local_label_type, network_population> labeled_population;
    for (const auto& s: pop) { labeled_population[s.label].emplace_back(s); }

    network_population unique_pop;

    for (auto& [_, pop]: labeled_population) {
        for (auto& a: pop) {
            if (a.begin < a.end)
                for (auto& b: pop) {
                    // check if intervals overlap and a != b
                    if ((a.begin != b.begin || a.end != b.end) &&
                        ((a.begin <= b.begin && a.end > b.begin) ||
                            (b.begin <= a.begin && b.end > a.begin))) {
                        // merge intervals and set one to empty
                        a.begin = std::min(a.begin, b.begin);
                        a.end = std::max(a.end, b.end);
                        b.begin = 0;
                        b.end = 0;
                    }
                }
        }
        for (const auto& a: pop) {
            // keep all non-empty entries
            if (a.begin < a.end) unique_pop.emplace_back(a);
        }
    }

    return unique_pop;
}

ARB_ARBOR_API network_population join(const network_population& a, network_population b) {
    b.insert(b.end(), a.begin(), a.end());
    return b;
}

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
        return uniform_rand_from_key_pair({unsigned(network_seed::selection), seed},
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
        return left.select(src_gid, src_label, dest_gid, dest_label) ||
               right.select(src_gid, src_label, dest_gid, dest_label);
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
        return !selection.select(src_gid, src_label, dest_gid, dest_label);
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

struct network_value::uniform_random_impl: public value_impl {
    unsigned seed = 0;
    std::array<double, 2> range;

    uniform_random_impl(unsigned rand_seed, std::array<double, 2> r): seed(rand_seed), range(r) {}

    double get(cell_gid_type src_gid,
        const cell_local_label_type& src_label,
        cell_gid_type dest_gid,
        const cell_local_label_type& dest_label) const override {
        if (range[0] > range[1]) return range[1];

        // random number between 0 and 1
        const auto rand_num = uniform_rand_from_key_pair({unsigned(network_seed::value), seed},
            hash_global_tag(src_gid, src_label.tag),
            hash_global_tag(dest_gid, dest_label.tag));

        return (range[1] - range[0]) * rand_num + range[0];
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

network_value network_value::uniform_random(unsigned seed, std::array<double, 2> range) {
    return {std::shared_ptr<value_impl>(new uniform_random_impl(seed, range))};
}

network_value network_value::custom(
    std::function<double(const cell_global_label_type&, const cell_global_label_type&)> func) {
    return {std::shared_ptr<value_impl>(new custom_impl(std::move(func)))};
}

network_value network_value::uniform(double value) {
    return {std::shared_ptr<value_impl>(new uniform_impl(value))};
}

cell_connection_network::cell_connection_network(network_value weight,
    network_value delay,
    network_selection selection,
    network_population src_pop,
    network_population dest_pop):
    selection_(std::move(selection)),
    weight_(std::move(weight)),
    delay_(std::move(delay)),
    src_pop_(std::move(src_pop)),
    dest_pop_(std::move(dest_pop)) {}

std::vector<cell_connection> cell_connection_network::generate(cell_gid_type gid) const {

    std::vector<cell_connection> connections;

    // For cell_connection, we only need to consider connections, where gid is
    // the destination
    for (const auto& dest: dest_pop_) {
        if (gid >= dest.begin && gid < dest.end) {
            cell_gid_type src_index = 0;
            for (const auto& src: src_pop_) {
                for (auto src_gid = src.begin; src_gid < src.end; ++src_gid, ++src_index) {
                    if (selection_.select(src_gid, src.label, gid, dest.label))
                        connections.emplace_back(cell_global_label_type{src_gid, src.label},
                            dest.label,
                            weight_.get(src_gid, src.label, gid, dest.label),
                            delay_.get(src_gid, src.label, gid, dest.label));
                }
            }
        }
    }

    // remove any duplicates
    std::sort(connections.begin(), connections.end());
    connections.erase(std::unique(connections.begin(), connections.end()), connections.end());

    return connections;
}

gap_junction_network::gap_junction_network(network_value weight,
    network_selection selection,
    network_population src_pop,
    network_population dest_pop):
    selection_(std::move(selection)),
    weight_(std::move(weight)),
    src_pop_(std::move(src_pop)),
    dest_pop_(std::move(dest_pop)) {}

std::vector<gap_junction_connection> gap_junction_network::generate(cell_gid_type gid) const {

    std::vector<gap_junction_connection> connections;

    // For gap junctions, we need to consider all connections, where gid is
    // the destination or the source
    for (const auto& dest: dest_pop_) {
        for (const auto& src: src_pop_) {
            for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                if (gid == dest_gid) {
                    for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                        if (selection_.select(src_gid, src.label, gid, dest.label))
                            connections.emplace_back(cell_global_label_type{src_gid, src.label},
                                dest.label,
                                weight_.get(src_gid, src.label, gid, dest.label));
                    }
                }
                if (gid >= src.begin && gid < src.end) {
                    if (selection_.select(dest_gid, dest.label, gid, src.label))
                        connections.emplace_back(cell_global_label_type{dest_gid, dest.label},
                            src.label,
                            weight_.get(dest_gid, dest.label, gid, src.label));
                }
            }
        }
    }

    // remove any duplicates
    std::sort(connections.begin(), connections.end());
    connections.erase(std::unique(connections.begin(), connections.end()), connections.end());

    return connections;
}

}  // namespace arb
