#pragma once

#include <arbor/common_types.hpp>
#include <arbor/math.hpp>

#include <vector>
#include <variant>
#include <optional>
#include <cstddef>
#include <type_traits>

namespace arb {

template <typename T, std::size_t DIM>
class spatial_tree {
public:
    static_assert(DIM >= 1, "Dimension of tree must be at least 1.");

    using value_type = T;
    using point_type = std::array<double, DIM>;
    using node_data = std::vector<spatial_tree>;
    using leaf_data = std::vector<std::pair<point_type, T>>;

    spatial_tree(): data_(leaf_data()) {}

    spatial_tree(leaf_data data): data_(std::move(data)) {}

    spatial_tree(std::size_t max_depth, std::size_t leaf_size_target, leaf_data data):
        data_(std::move(data)) {
        if (std::get<leaf_data>(data_).size()) {
            min_.fill(std::numeric_limits<double>::max());
            max_.fill(std::numeric_limits<double>::min());

            for (const auto &[p, _]: std::get<leaf_data>(data_)) {
                for (std::size_t i = 0; i < DIM; ++i) {
                    if (p[i] < min_[i]) min_[i] = p[i];
                    if (p[i] > max_[i]) max_[i] = p[i];
                }
            }

            for (std::size_t i = 0; i < DIM; ++i) { mid_[i] = (max_[i] - min_[i]) / 2.0 + min_[i]; }

            split_recursively(0, max_depth, leaf_size_target);
        }
    }

    // Iterate over all points recursively.
    // func must have signature `void func(const point_type&, const T&)`.
    template <typename F>
    inline void for_each(const F &func) const {
        std::visit(
            [&](auto &&arg) {
                using arg_type = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<arg_type, node_data>) {
                    for (const auto &node: arg) { node.for_each(func); }
                }
                if constexpr (std::is_same_v<arg_type, leaf_data>) {
                    for (const auto &[p, d]: arg) { func(p, d); }
                }
            },
            data_);
    }

    // Iterate over all points within the given bounding box recursively.
    // func must have signature `void func(const point_type&, const T&)`.
    template <typename F>
    inline void bounding_box_for_each(const point_type &box_min,
        const point_type &box_max,
        const F &func) const {
        auto all_smaller_eq = [](const point_type &lhs, const point_type &rhs) {
            bool result = true;
            for (std::size_t i = 0; i < DIM; ++i) { result &= lhs[i] <= rhs[i]; }
            return result;
        };

        std::visit(
            [&](auto &&arg) {
                using arg_type = std::decay_t<decltype(arg)>;

                if (all_smaller_eq(box_min, min_) && all_smaller_eq(max_, box_max)) {
                    // sub-nodes fully inside box -> call without further boundary
                    // checks

                    if constexpr (std::is_same_v<arg_type, node_data>) {
                        for (const auto &node: arg) { node.template for_each<F>(func); }
                    }
                    if constexpr (std::is_same_v<arg_type, leaf_data>) {
                        for (const auto &[p, d]: arg) { func(p, d); }
                    }
                }
                else {
                    // sub-nodes partially overlap bounding box
                    if constexpr (std::is_same_v<arg_type, node_data>) {
                        for (const auto &node: arg) {
                            if (all_smaller_eq(node.min_, box_max) &&
                                all_smaller_eq(box_min, node.max_))
                                node.template bounding_box_for_each<F>(box_min, box_max, func);
                        }
                    }
                    if constexpr (std::is_same_v<arg_type, leaf_data>) {
                        for (const auto &[p, d]: arg) {
                            if (all_smaller_eq(p, box_max) && all_smaller_eq(box_min, p)) {
                                func(p, d);
                            }
                        }
                    }
                }
            },
            data_);
    }

private:
    // The index of the sub node containing p
    inline std::size_t sub_node_index(const point_type &p) {
        std::size_t index = 0;
        for (std::size_t i = 0; i < DIM; ++i) { index += i * 2 * (p[i] >= mid_[i]); }
        return index;
    }

    // converts a leaf node into a tree node by splitting the data_ into 8
    // sub-nodes
    inline void split_current_node() {
        if (std::holds_alternative<leaf_data>(data_)) {
            constexpr auto divisor = math::pow<std::size_t, std::size_t>(2, DIM);

            node_data new_nodes;
            new_nodes.reserve(divisor);

            // assign each point to sub-node
            std::array<leaf_data, divisor> new_leaf_data;
            for (const auto &[p, d]: std::get<leaf_data>(data_)) {
                new_leaf_data[sub_node_index(p)].emplace_back(p, d);
            }

            // move sorted data_ into new sub-nodes
            for (auto &d: new_leaf_data) { new_nodes.emplace_back(std::move(d)); }

            // replace current data_ with new sub-nodes
            this->data_ = std::move(new_nodes);
        }
    }

    inline void split_recursively(std::size_t depth,
        std::size_t max_depth,
        std::size_t leaf_size_target) {
        if (depth >= max_depth) {
            // if maximum depth is reached do nothing for leaf nodes or combine all
            // sub-node data for non-leaf nodes
            std::visit(
                [&](auto &&arg) {
                    using arg_type = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<arg_type, node_data>) {
                        leaf_data combined_data;
                        this->for_each([&](const point_type &p, const T &d) {
                            combined_data.emplace_back(p, d);
                        });
                        this->data_ = std::move(combined_data);
                    }
                },
                data_);
        }
        else {

            std::visit(
                [&](auto &&arg) {
                    using arg_type = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<arg_type, node_data>) {
                        for (auto &node: arg) {
                            node.split_recursively(depth + 1, max_depth, leaf_size_target);
                        }
                    }
                    if constexpr (std::is_same_v<arg_type, leaf_data>) {
                        if (arg.size() > leaf_size_target) {
                            // converts from leaf to node. Invalidates arg.
                            this->split_current_node();
                            // call recursive split again for same depth, which will execute
                            // the other branch
                            this->split_recursively(depth, max_depth, leaf_size_target);
                        }
                    }
                },
                data_);
        }
    }

    point_type min_, max_, mid_;
    std::variant<node_data, leaf_data> data_;
};

}  // namespace arb
