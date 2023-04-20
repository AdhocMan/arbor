#include <arbor/common_types.hpp>
#include <arbor/network.hpp>

#include <Random123/threefry.h>
#include <Random123/boxmuller.hpp>
#include <Random123/uniform.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <type_traits>
#include <variant>
#include <vector>

#include "network_impl.hpp"

namespace arb {

namespace {

// Partial seed to use for network_value and network_selection generation.
// Different seed for each type to avoid unintentional correlation.
enum class network_seed : unsigned {
    selection_random = 2058443,
    value_uniform = 48202,
    value_normal = 8405,
    value_truncated_normal = 380237,
    site_info = 984293
};

// We only need minimal hash collisions and good spread over the hash range, because this will be
// used as input for random123, which then provides all desired hash properties.
// std::hash is implementation dependent, so we define our own for reproducibility.

std::uint64_t simple_string_hash(const std::string_view& s) {
    // use fnv1a hash algorithm
    constexpr std::uint64_t prime = 1099511628211ull;
    std::uint64_t h = 14695981039346656037ull;

    for (auto c: s) {
        h ^= c;
        h *= prime;
    }

    return h;
}

double uniform_rand_from_key_pair(std::array<unsigned, 2> seed,
    network_hash_type key_a,
    network_hash_type key_b) {
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

struct network_selection_all_impl: public network_selection_impl {
    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return true;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override { os << "(all)"; }
};

struct network_selection_none_impl: public network_selection_impl {

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return false;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return false;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return false;
    }

    void print(std::ostream& os) const override { os << "(none)"; }
};

struct network_selection_source_cell_kind_impl: public network_selection_impl {
    cell_kind select_kind;

    explicit network_selection_source_cell_kind_impl(cell_kind k): select_kind(k) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return src.kind == select_kind;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return kind == select_kind;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override {
        os << "(source-cell-kind (";
        switch (select_kind) {
        case arb::cell_kind::spike_source: os << "spike-source"; break;
        case arb::cell_kind::cable: os << "cable"; break;
        case arb::cell_kind::lif: os << "lif"; break;
        case arb::cell_kind::benchmark: os << "benchmark"; break;
        }
        os << "-cell))";
    }
};

struct network_selection_destination_cell_kind_impl: public network_selection_impl {
    cell_kind select_kind;

    explicit network_selection_destination_cell_kind_impl(cell_kind k): select_kind(k) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return dest.kind == select_kind;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return kind == select_kind;
    }

    void print(std::ostream& os) const override {
        os << "(destination-cell-kind (";
        switch (select_kind) {
        case arb::cell_kind::spike_source: os << "spike-source"; break;
        case arb::cell_kind::cable: os << "cable"; break;
        case arb::cell_kind::lif: os << "lif"; break;
        case arb::cell_kind::benchmark: os << "benchmark"; break;
        }
        os << "-cell))";
    }
};

struct network_selection_source_label_impl: public network_selection_impl {
    std::vector<cell_tag_type> sorted_labels;

    explicit network_selection_source_label_impl(std::vector<cell_tag_type> labels):
        sorted_labels(std::move(labels)) {
        std::sort(sorted_labels.begin(), sorted_labels.end());
    }

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return std::binary_search(sorted_labels.begin(), sorted_labels.end(), src.label);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return std::binary_search(sorted_labels.begin(), sorted_labels.end(), label);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override {
        os << "(source-label";
        for (const auto& l: sorted_labels) { os << " \"" << l << "\""; }
        os << ")";
    }
};

struct network_selection_destination_label_impl: public network_selection_impl {
    std::vector<cell_tag_type> sorted_labels;

    explicit network_selection_destination_label_impl(std::vector<cell_tag_type> labels):
        sorted_labels(std::move(labels)) {
        std::sort(sorted_labels.begin(), sorted_labels.end());
    }

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return std::binary_search(sorted_labels.begin(), sorted_labels.end(), dest.label);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return std::binary_search(sorted_labels.begin(), sorted_labels.end(), label);
    }

    void print(std::ostream& os) const override {
        os << "(destination-label";
        for (const auto& l: sorted_labels) { os << " \"" << l << "\""; }
        os << ")";
    }
};

struct network_selection_source_cell_impl: public network_selection_impl {
    std::vector<cell_gid_type> sorted_gids;

    explicit network_selection_source_cell_impl(std::vector<cell_gid_type> gids):
        sorted_gids(std::move(gids)) {
        std::sort(sorted_gids.begin(), sorted_gids.end());
    }

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return std::binary_search(sorted_gids.begin(), sorted_gids.end(), src.gid);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return std::binary_search(sorted_gids.begin(), sorted_gids.end(), gid);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override {
        os << "(source-cell";
        for (const auto& g: sorted_gids) { os << " " << g; }
        os << ")";
    }
};

struct network_selection_source_cell_range_impl: public network_selection_impl {
    cell_gid_type gid_begin, gid_end, step;

    network_selection_source_cell_range_impl(gid_range r):
        gid_begin(r.begin),
        gid_end(r.end),
        step(r.step) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return src.gid >= gid_begin && src.gid < gid_end && !((src.gid - gid_begin) % step);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return gid >= gid_begin && gid < gid_end && !((gid - gid_begin) % step);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override {
        os << "(source-cell (gid-range " << gid_begin << " " << gid_end << " " << step << "))";
    }
};

struct network_selection_destination_cell_impl: public network_selection_impl {
    std::vector<cell_gid_type> sorted_gids;

    network_selection_destination_cell_impl(std::vector<cell_gid_type> gids):
        sorted_gids(std::move(gids)) {
        std::sort(sorted_gids.begin(), sorted_gids.end());
    }

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return std::binary_search(sorted_gids.begin(), sorted_gids.end(), dest.gid);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return std::binary_search(sorted_gids.begin(), sorted_gids.end(), gid);
    }

    void print(std::ostream& os) const override {
        os << "(destination-cell";
        for (const auto& g: sorted_gids) { os << " " << g; }
        os << ")";
    }
};

struct network_selection_destination_cell_range_impl: public network_selection_impl {
    cell_gid_type gid_begin, gid_end, step;

    network_selection_destination_cell_range_impl(gid_range r):
        gid_begin(r.begin),
        gid_end(r.end),
        step(r.step) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return dest.gid >= gid_begin && dest.gid < gid_end && !((dest.gid - gid_begin) % step);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return gid >= gid_begin && gid < gid_end && !((gid - gid_begin) % step);
    }

    void print(std::ostream& os) const override {
        os << "(destination-cell (gid-range " << gid_begin << " " << gid_end << " " << step << "))";
    }
};

struct network_selection_chain_impl: public network_selection_impl {
    std::vector<cell_gid_type> gids;  // preserved order of ring
    std::vector<cell_gid_type> sorted_gids;
    network_selection_chain_impl(std::vector<cell_gid_type> gids): gids(std::move(gids)) {
        sorted_gids = this->gids;  // copy
        std::sort(sorted_gids.begin(), sorted_gids.end());
    }

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        if (gids.empty()) return false;

        // gids size always > 0 frome here on

        // First check if both are part of ring
        if (!std::binary_search(sorted_gids.begin(), sorted_gids.end(), src.gid) ||
            !std::binary_search(sorted_gids.begin(), sorted_gids.end(), dest.gid))
            return false;

        for (std::size_t i = 0; i < gids.size() - 1; ++i) {
            // return true if neighbors in gids list
            if ((src.gid == gids[i] && dest.gid == gids[i + 1])) return true;
        }

        return false;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return !sorted_gids.empty() &&
               std::binary_search(sorted_gids.begin(), sorted_gids.end() - 1, gid);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return !sorted_gids.empty() &&
               std::binary_search(sorted_gids.begin() + 1, sorted_gids.end(), gid);
    }

    void print(std::ostream& os) const override {
        os << "(chain";
        for (const auto& g: gids) { os << " " << g; }
        os << ")";
    }
};

struct network_selection_chain_range_impl: public network_selection_impl {
    cell_gid_type gid_begin, gid_end, step;

    network_selection_chain_range_impl(gid_range r):
        gid_begin(r.begin),
        gid_end(r.end),
        step(r.step) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        if (src.gid < gid_begin || src.gid >= gid_end || dest.gid < gid_begin ||
            dest.gid >= gid_end)
            return false;

        return src.gid + step == dest.gid && !((src.gid - gid_begin) % step);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        // Return false if outside range or if equal to last element, which cannot be a source
        if (gid < gid_begin || gid >= gid_end - 1) return false;
        return !((gid - gid_begin) % step);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        // Return false if outside range or if equal to first element, which cannot be a destination
        if (gid <= gid_begin || gid >= gid_end) return false;
        return !((gid - gid_begin) % step);
    }

    void print(std::ostream& os) const override {
        os << "(chain (gid-range " << gid_begin << " " << gid_end << " " << step << "))";
    }
};

struct network_selection_reverse_chain_range_impl: public network_selection_impl {
    cell_gid_type gid_begin, gid_end, step;

    network_selection_reverse_chain_range_impl(gid_range r):
        gid_begin(r.begin),
        gid_end(r.end),
        step(r.step) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        if (src.gid < gid_begin || src.gid >= gid_end || dest.gid < gid_begin ||
            dest.gid >= gid_end)
            return false;

        return dest.gid + step == src.gid && !((src.gid - gid_begin) % step);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        // Return false if outside range or if equal to first element, which cannot be a source
        if (gid <= gid_begin || gid >= gid_end) return false;
        return !((gid - gid_begin) % step);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        // Return false if outside range or if equal to last element, which cannot be a destination
        if (gid < gid_begin || gid >= gid_end - 1) return false;
        return !((gid - gid_begin) % step);
    }

    void print(std::ostream& os) const override {
        os << "(chain-reverse (gid-range " << gid_begin << " " << gid_end << " " << step << "))";
    }
};

struct network_selection_complement_impl: public network_selection_impl {
    std::shared_ptr<network_selection_impl> selection;

    explicit network_selection_complement_impl(std::shared_ptr<network_selection_impl> s):
        selection(std::move(s)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return !selection->select_connection(src, dest);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;  // cannot exclude any because source selection cannot be complemented without
                      // knowing selection criteria.
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;  // cannot exclude any because destination selection cannot be complemented
                      // without knowing selection criteria.
    }

    void initialize(const network_label_dict& dict) override { selection->initialize(dict); };

    void print(std::ostream& os) const override {
        os << "(complement ";
        selection->print(os);
        os << ")";
    }
};

struct network_selection_named_impl: public network_selection_impl {
    using impl_pointer_type = std::shared_ptr<network_selection_impl>;

    impl_pointer_type selection;
    std::string selection_name;

    explicit network_selection_named_impl(std::string name): selection_name(std::move(name)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        if (!selection)
            throw arbor_internal_error("Trying to use unitialized named network selection.");
        return selection->select_connection(src, dest);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        if (!selection)
            throw arbor_internal_error("Trying to use unitialized named network selection.");
        return selection->select_source(kind, gid, label);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        if (!selection)
            throw arbor_internal_error("Trying to use unitialized named network selection.");
        return selection->select_destination(kind, gid, label);
    }

    void initialize(const network_label_dict& dict) override {
        auto s = dict.selection(selection_name);
        if (!s.has_value())
            throw arbor_exception(
                std::string("Network selection with label \"") + selection_name + "\" not found.");
        selection = thingify(s.value(), dict);
    };

    void print(std::ostream& os) const override {
        os << "(network-selection \"" << selection_name << "\")";
    }
};

struct network_selection_inter_cell_impl: public network_selection_impl {
    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return src.gid != dest.gid;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override { os << "(inter-cell)"; }
};

struct network_selection_custom_impl: public network_selection_impl {
    network_selection::custom_func_type func;

    explicit network_selection_custom_impl(network_selection::custom_func_type f):
        func(std::move(f)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return func(src, dest);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override { os << "(custom-network-selection)"; }
};

struct network_selection_distance_lt_impl: public network_selection_impl {
    double d;

    explicit network_selection_distance_lt_impl(double d): d(d) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return distance(src.global_location, dest.global_location) < d;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    std::optional<double> max_distance() const override { return d; }

    void print(std::ostream& os) const override { os << "(distance-lt " << d << ")"; }
};

struct network_selection_distance_gt_impl: public network_selection_impl {
    double d;

    explicit network_selection_distance_gt_impl(double d): d(d) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return distance(src.global_location, dest.global_location) > d;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void print(std::ostream& os) const override { os << "(distance-gt " << d << ")"; }
};

struct network_selection_random_impl: public network_selection_impl {
    unsigned seed;

    network_value p_value;
    std::shared_ptr<network_value_impl> probability; // may be null if unitialize(...) not called

    network_selection_random_impl(unsigned seed, network_value p): seed(seed), p_value(std::move(p)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        if (!probability)
            throw arbor_internal_error("Trying to use unitialized named network selection.");
        const auto r = uniform_rand_from_key_pair(
            {unsigned(network_seed::selection_random), seed}, src.hash, dest.hash);
        const auto p = (probability->get(src, dest));
        return r < p;
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return true;
    }

    void initialize(const network_label_dict& dict) override {
        probability = thingify(p_value, dict);
    };

    void print(std::ostream& os) const override {
        os << "(random " << seed << " ";
        os << p_value;
        os << ")";
    }
};

struct network_selection_intersect_impl: public network_selection_impl {
    std::shared_ptr<network_selection_impl> left, right;

    network_selection_intersect_impl(std::shared_ptr<network_selection_impl> l,
        std::shared_ptr<network_selection_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return left->select_connection(src, dest) && right->select_connection(src, dest);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_source(kind, gid, label) && right->select_source(kind, gid, label);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_destination(kind, gid, label) &&
               right->select_destination(kind, gid, label);
    }

    std::optional<double> max_distance() const override {
        const auto d_left = left->max_distance();
        const auto d_right = right->max_distance();

        if (d_left && d_right) return std::min(d_left.value(), d_right.value());
        if (d_left) return d_left.value();
        if (d_right) return d_right.value();

        return std::nullopt;
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(intersect ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_selection_join_impl: public network_selection_impl {
    std::shared_ptr<network_selection_impl> left, right;

    network_selection_join_impl(std::shared_ptr<network_selection_impl> l,
        std::shared_ptr<network_selection_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return left->select_connection(src, dest) || right->select_connection(src, dest);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_source(kind, gid, label) || right->select_source(kind, gid, label);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_destination(kind, gid, label) ||
               right->select_destination(kind, gid, label);
    }

    std::optional<double> max_distance() const override {
        const auto d_left = left->max_distance();
        const auto d_right = right->max_distance();

        if (d_left && d_right) return std::max(d_left.value(), d_right.value());

        return std::nullopt;
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(join ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_selection_symmetric_difference_impl: public network_selection_impl {
    std::shared_ptr<network_selection_impl> left, right;

    network_selection_symmetric_difference_impl(std::shared_ptr<network_selection_impl> l,
        std::shared_ptr<network_selection_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return left->select_connection(src, dest) ^ right->select_connection(src, dest);
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_source(kind, gid, label) || right->select_source(kind, gid, label);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_destination(kind, gid, label) ||
               right->select_destination(kind, gid, label);
    }

    std::optional<double> max_distance() const override {
        const auto d_left = left->max_distance();
        const auto d_right = right->max_distance();

        if (d_left && d_right) return std::max(d_left.value(), d_right.value());

        return std::nullopt;
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(symmetric-difference ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_selection_difference_impl: public network_selection_impl {
    std::shared_ptr<network_selection_impl> left, right;

    network_selection_difference_impl(std::shared_ptr<network_selection_impl> l,
        std::shared_ptr<network_selection_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    bool select_connection(const network_site_info& src,
        const network_site_info& dest) const override {
        return left->select_connection(src, dest) && !(right->select_connection(src, dest));
    }

    bool select_source(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_source(kind, gid, label);
    }

    bool select_destination(cell_kind kind,
        cell_gid_type gid,
        const std::string_view& label) const override {
        return left->select_destination(kind, gid, label);
    }

    std::optional<double> max_distance() const override {
        const auto d_left = left->max_distance();

        if (d_left) return d_left.value();

        return std::nullopt;
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(difference ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_value_scalar_impl: public network_value_impl {
    double value;

    network_value_scalar_impl(double v): value(v) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return value;
    }

    void print(std::ostream& os) const override { os << "(scalar " << value << ")"; }
};


struct network_value_distance_impl: public network_value_impl {
    double scale;

    network_value_distance_impl(double s): scale(s) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return scale * distance(src.global_location, dest.global_location);
    }

    void print(std::ostream& os) const override { os << "(distance " << scale << ")"; }
};

struct network_value_uniform_distribution_impl: public network_value_impl {
    unsigned seed = 0;
    std::array<double, 2> range;

    network_value_uniform_distribution_impl(unsigned rand_seed, const std::array<double, 2>& r):
        seed(rand_seed),
        range(r) {
        if (range[0] >= range[1])
            throw std::invalid_argument("Uniform distribution: invalid range");
    }

    double get(const network_site_info& src, const network_site_info& dest) const override {
        if (range[0] > range[1]) return range[1];

        // random number between 0 and 1
        const auto rand_num = uniform_rand_from_key_pair(
            {unsigned(network_seed::value_uniform), seed}, src.hash, dest.hash);

        return (range[1] - range[0]) * rand_num + range[0];
    }

    void print(std::ostream& os) const override {
        os << "(uniform-distribution " << seed << " " << range[0] << " " << range[1] << ")";
    }
};

struct network_value_normal_distribution_impl: public network_value_impl {
    unsigned seed = 0;
    double mean = 0.0;
    double std_deviation = 1.0;

    network_value_normal_distribution_impl(unsigned rand_seed, double mean_, double std_deviation_):
        seed(rand_seed),
        mean(mean_),
        std_deviation(std_deviation_) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return mean + std_deviation *
                          normal_rand_from_key_pair(
                              {unsigned(network_seed::value_normal), seed}, src.hash, dest.hash);
    }

    void print(std::ostream& os) const override {
        os << "(normal-distribution " << seed << " " << mean << " " << std_deviation << ")";
    }
};

struct network_value_truncated_normal_distribution_impl: public network_value_impl {
    unsigned seed = 0;
    double mean = 0.0;
    double std_deviation = 1.0;
    std::array<double, 2> range;

    network_value_truncated_normal_distribution_impl(unsigned rand_seed,
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

    double get(const network_site_info& src, const network_site_info& dest) const override {

        const auto src_hash = src.hash;
        auto dest_hash = dest.hash;

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

    void print(std::ostream& os) const override {
        os << "(truncated-normal-distribution " << seed << " " << mean << " " << std_deviation
           << " " << range[0] << " " << range[1] << ")";
    }
};

struct network_value_custom_impl: public network_value_impl {
    network_value::custom_func_type func;

    network_value_custom_impl(network_value::custom_func_type f): func(std::move(f)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return func(src, dest);
    }

    void print(std::ostream& os) const override { os << "(custom-network-value)"; }
};

struct network_value_named_impl: public network_value_impl {
    using impl_pointer_type = std::shared_ptr<network_value_impl>;

    impl_pointer_type value;
    std::string value_name;

    explicit network_value_named_impl(std::string name): value_name(std::move(name)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        if (!value) throw arbor_internal_error("Trying to use unitialized named network value.");
        return value->get(src, dest);
    }

    void initialize(const network_label_dict& dict) override {
        auto v = dict.value(value_name);
        if (!v.has_value())
            throw arbor_exception(
                std::string("Network value with label \"") + value_name + "\" not found.");
        value = thingify(v.value(), dict);
    };

    void print(std::ostream& os) const override {
        os << "(network-value \"" << value_name << "\")";
    }
};

struct network_value_add_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> left, right;

    network_value_add_impl(std::shared_ptr<network_value_impl> l,
        std::shared_ptr<network_value_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return left->get(src, dest) + right->get(src, dest);
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(add ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};


struct network_value_mul_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> left, right;

    network_value_mul_impl(std::shared_ptr<network_value_impl> l,
        std::shared_ptr<network_value_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return left->get(src, dest) * right->get(src, dest);
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(mul ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};


struct network_value_sub_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> left, right;

    network_value_sub_impl(std::shared_ptr<network_value_impl> l,
        std::shared_ptr<network_value_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return left->get(src, dest) - right->get(src, dest);
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(sub ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_value_div_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> left, right;

    network_value_div_impl(std::shared_ptr<network_value_impl> l,
        std::shared_ptr<network_value_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        const auto v_right = right ->get(src,dest);
        if (!v_right) throw arbor_exception("network_value: division by 0.");
        return left->get(src, dest) / right->get(src, dest);
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(div ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_value_max_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> left, right;

    network_value_max_impl(std::shared_ptr<network_value_impl> l,
        std::shared_ptr<network_value_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return std::max(left->get(src, dest), right->get(src, dest));
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(max ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_value_min_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> left, right;

    network_value_min_impl(std::shared_ptr<network_value_impl> l,
        std::shared_ptr<network_value_impl> r):
        left(std::move(l)),
        right(std::move(r)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return std::min(left->get(src, dest), right->get(src, dest));
    }

    void initialize(const network_label_dict& dict) override {
        left->initialize(dict);
        right->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(min ";
        left->print(os);
        os << " ";
        right->print(os);
        os << ")";
    }
};

struct network_value_exp_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> value;

    network_value_exp_impl(std::shared_ptr<network_value_impl> v):
        value(std::move(v)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        return std::exp(value->get(src, dest));
    }

    void initialize(const network_label_dict& dict) override {
        value->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(exp ";
        value->print(os);
        os << ")";
    }
};

struct network_value_log_impl: public network_value_impl {
    std::shared_ptr<network_value_impl> value;

    network_value_log_impl(std::shared_ptr<network_value_impl> v):
        value(std::move(v)) {}

    double get(const network_site_info& src, const network_site_info& dest) const override {
        const auto v = value->get(src, dest);
        if (v <= 0.0) throw arbor_exception("network_value: log of value <= 0.0.");
        return std::log(value->get(src, dest));
    }

    void initialize(const network_label_dict& dict) override {
        value->initialize(dict);
    };

    void print(std::ostream& os) const override {
        os << "(log ";
        value->print(os);
        os << ")";
    }
};

}  // namespace

network_site_info::network_site_info(cell_gid_type gid,
    cell_lid_type lid,
    cell_kind kind,
    std::string_view label,
    mlocation location,
    mpoint global_location):
    gid(gid),
    lid(lid),
    kind(kind),
    label(std::move(label)),
    location(location),
    global_location(global_location) {

    std::uint64_t label_hash = simple_string_hash(this->label);
    static_assert(sizeof(decltype(mlocation::pos)) == sizeof(std::uint64_t));
    std::uint64_t loc_pos_hash = *reinterpret_cast<const std::uint64_t*>(&location.pos);

    const auto seed = static_cast<std::uint64_t>(network_seed::site_info);

    using rand_type = r123::Threefry4x64;
    const rand_type::ctr_type seed_input = {{seed, 2 * seed, 3 * seed, 4 * seed}};
    const rand_type::key_type key = {{gid, label_hash, location.branch, loc_pos_hash}};

    rand_type gen;
    hash = gen(seed_input, key)[0];
}

network_selection::network_selection(std::shared_ptr<network_selection_impl> impl):
    impl_(std::move(impl)) {}

network_selection network_selection::intersect(network_selection left, network_selection right) {
    return network_selection(std::make_shared<network_selection_intersect_impl>(
        std::move(left.impl_), std::move(right.impl_)));
}

network_selection network_selection::join(network_selection left, network_selection right) {
    return network_selection(std::make_shared<network_selection_join_impl>(
        std::move(left.impl_), std::move(right.impl_)));
}

network_selection network_selection::symmetric_difference(network_selection left,
    network_selection right) {
    return network_selection(std::make_shared<network_selection_symmetric_difference_impl>(
        std::move(left.impl_), std::move(right.impl_)));
}

network_selection network_selection::difference(network_selection left, network_selection right) {
    return network_selection(std::make_shared<network_selection_difference_impl>(
        std::move(left.impl_), std::move(right.impl_)));
}

network_selection network_selection::all() {
    return network_selection(std::make_shared<network_selection_all_impl>());
}

network_selection network_selection::none() {
    return network_selection(std::make_shared<network_selection_none_impl>());
}

network_selection network_selection::named(std::string name) {
    return network_selection(std::make_shared<network_selection_named_impl>(std::move(name)));
}

network_selection network_selection::source_cell_kind(cell_kind kind) {
    return network_selection(std::make_shared<network_selection_source_cell_kind_impl>(kind));
}

network_selection network_selection::destination_cell_kind(cell_kind kind) {
    return network_selection(std::make_shared<network_selection_destination_cell_kind_impl>(kind));
}

network_selection network_selection::source_label(std::vector<cell_tag_type> labels) {
    return network_selection(
        std::make_shared<network_selection_source_label_impl>(std::move(labels)));
}

network_selection network_selection::destination_label(std::vector<cell_tag_type> labels) {
    return network_selection(
        std::make_shared<network_selection_destination_label_impl>(std::move(labels)));
}

network_selection network_selection::source_cell(std::vector<cell_gid_type> gids) {
    return network_selection(std::make_shared<network_selection_source_cell_impl>(std::move(gids)));
}

network_selection network_selection::source_cell(gid_range range) {
    return network_selection(std::make_shared<network_selection_source_cell_range_impl>(range));
}

network_selection network_selection::destination_cell(std::vector<cell_gid_type> gids) {
    return network_selection(
        std::make_shared<network_selection_destination_cell_impl>(std::move(gids)));
}

network_selection network_selection::destination_cell(gid_range range) {
    return network_selection(
        std::make_shared<network_selection_destination_cell_range_impl>(range));
}

network_selection network_selection::chain(std::vector<cell_gid_type> gids) {
    return network_selection(std::make_shared<network_selection_chain_impl>(std::move(gids)));
}

network_selection network_selection::chain(gid_range range) {
    return network_selection(std::make_shared<network_selection_chain_range_impl>(range));
}

network_selection network_selection::chain_reverse(gid_range range) {
    return network_selection(std::make_shared<network_selection_reverse_chain_range_impl>(range));
}

network_selection network_selection::complement(network_selection s) {
    return network_selection(
        std::make_shared<network_selection_complement_impl>(std::move(s.impl_)));
}

network_selection network_selection::inter_cell() {
    return network_selection(std::make_shared<network_selection_inter_cell_impl>());
}

network_selection network_selection::random(unsigned seed, network_value p) {
    return network_selection(
        std::make_shared<network_selection_random_impl>(seed, std::move(p)));
}

network_selection network_selection::custom(custom_func_type func) {
    return network_selection(std::make_shared<network_selection_custom_impl>(std::move(func)));
}

network_selection network_selection::distance_lt(double d) {
    return network_selection(std::make_shared<network_selection_distance_lt_impl>(d));
}

network_selection network_selection::distance_gt(double d) {
    return network_selection(std::make_shared<network_selection_distance_gt_impl>(d));
}

network_value::network_value(std::shared_ptr<network_value_impl> impl): impl_(std::move(impl)) {}

network_value network_value::scalar(double value) {
    return network_value(std::make_shared<network_value_scalar_impl>(value));
}

network_value network_value::distance(double scale) {
    return network_value(std::make_shared<network_value_distance_impl>(scale));
}

network_value network_value::uniform_distribution(unsigned seed,
    const std::array<double, 2>& range) {
    return network_value(std::make_shared<network_value_uniform_distribution_impl>(seed, range));
}

network_value network_value::normal_distribution(unsigned seed, double mean, double std_deviation) {
    return network_value(
        std::make_shared<network_value_normal_distribution_impl>(seed, mean, std_deviation));
}

network_value network_value::truncated_normal_distribution(unsigned seed,
    double mean,
    double std_deviation,
    const std::array<double, 2>& range) {
    return network_value(std::make_shared<network_value_truncated_normal_distribution_impl>(
        seed, mean, std_deviation, range));
}

network_value network_value::custom(custom_func_type func) {
    return network_value(std::make_shared<network_value_custom_impl>(std::move(func)));
}

network_value network_value::named(std::string name) {
    return network_value(std::make_shared<network_value_named_impl>(std::move(name)));
}

network_label_dict& network_label_dict::set(const std::string& name, network_selection s) {
    selections_.insert_or_assign(name, std::move(s));
    return *this;
}

network_label_dict& network_label_dict::set(const std::string& name, network_value v) {
    values_.insert_or_assign(name, std::move(v));
    return *this;
}

network_value network_value::add(network_value left, network_value right) {
    return network_value(
        std::make_shared<network_value_add_impl>(std::move(left.impl_), std::move(right.impl_)));
}

network_value network_value::sub(network_value left, network_value right) {
    return network_value(
        std::make_shared<network_value_sub_impl>(std::move(left.impl_), std::move(right.impl_)));
}

network_value network_value::mul(network_value left, network_value right) {
    return network_value(
        std::make_shared<network_value_mul_impl>(std::move(left.impl_), std::move(right.impl_)));
}

network_value network_value::div(network_value left, network_value right) {
    return network_value(
        std::make_shared<network_value_div_impl>(std::move(left.impl_), std::move(right.impl_)));
}

network_value network_value::exp(network_value v) {
    return network_value(std::make_shared<network_value_exp_impl>(std::move(v.impl_)));
}

network_value network_value::log(network_value v) {
    return network_value(std::make_shared<network_value_log_impl>(std::move(v.impl_)));
}

network_value network_value::min(network_value left, network_value right) {
    return network_value(
        std::make_shared<network_value_min_impl>(std::move(left.impl_), std::move(right.impl_)));
}

network_value network_value::max(network_value left, network_value right) {
    return network_value(
        std::make_shared<network_value_max_impl>(std::move(left.impl_), std::move(right.impl_)));
}

std::optional<network_selection> network_label_dict::selection(const std::string& name) const {
    auto it = selections_.find(name);
    if (it != selections_.end()) return it->second;

    return std::nullopt;
}

std::optional<network_value> network_label_dict::value(const std::string& name) const {
    auto it = values_.find(name);
    if (it != values_.end()) return it->second;

    return std::nullopt;
}

ARB_ARBOR_API std::ostream& operator<<(std::ostream& os, const network_selection& s) {
    if (s.impl_) s.impl_->print(os);
    return os;
}

ARB_ARBOR_API std::ostream& operator<<(std::ostream& os, const network_value& v) {
    if (v.impl_) v.impl_->print(os);
    return os;
}

}  // namespace arb
