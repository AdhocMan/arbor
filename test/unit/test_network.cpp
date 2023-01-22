#include <gtest/gtest.h>

#include <arbor/network.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <unordered_map>

using namespace arb;

template <typename Group_t>
static bool is_gj_population_member(const std::vector<Group_t>& pop,
    const cell_global_label_type& l) {
    for (const auto& range: pop) {
        if (range.gid_begin <= l.gid && range.gid_end > l.gid) {
            for (const auto& group_label: range.labels) {
                if (group_label == l.label) return true;
            }
        }
    }
    return false;
}

template <typename Group_t>
static bool is_src_population_member(const std::vector<Group_t>& pop,
    const cell_global_label_type& l) {
    for (const auto& range: pop) {
        if (range.gid_begin <= l.gid && range.gid_end > l.gid) {
            for (const auto& group_label: range.src_labels) {
                if (group_label == l.label) return true;
            }
        }
    }
    return false;
}
template <typename Group_t>
static bool is_dest_population_member(const std::vector<Group_t>& pop,
    const cell_global_label_type& l) {
    for (const auto& range: pop) {
        if (range.gid_begin <= l.gid && range.gid_end > l.gid) {
            for (const auto& group_label: range.dest_labels) {
                if (group_label == l.label) return true;
            }
        }
    }
    return false;
}

template <typename Weight_t, typename Delay_t, typename Select_t, typename Group_t>
static long check_cell_connections(Weight_t weight,
    Delay_t delay,
    Select_t selector,
    const std::vector<Group_t>& pop) {

    cell_connection_network cn(weight, delay, selector, pop);

    // A gid may occur more than once, so gather unique sizes first
    std::unordered_map<cell_gid_type, std::size_t> sizes;

    for (const auto& group: pop) {
        for (auto gid = group.gid_begin; gid < group.gid_end; ++gid) {
            auto connections = cn.generate(gid);
            sizes.insert_or_assign(gid, connections.size());
            for (const auto& c: connections) {
                EXPECT_TRUE(is_src_population_member(pop, c.source));
                EXPECT_TRUE(is_dest_population_member(pop, {gid, c.dest}));
                EXPECT_TRUE(selector(c.source, {gid, c.dest}));
                EXPECT_EQ(weight(c.source, {gid, c.dest}), c.weight);
                EXPECT_EQ(delay(c.source, {gid, c.dest}), c.delay);
            }
        }
    }

    return std::accumulate(
        sizes.begin(), sizes.end(), 0, [](const auto& a, const auto& b) { return a + b.second; });
}

template <typename Weight_t, typename Select_t, typename Group_t>
static long check_gap_junctions(Weight_t weight,
    Select_t selector,
    const std::vector<Group_t>& pop) {

    gap_junction_network gn(weight, selector, pop);

    std::unordered_map<cell_gid_type, std::vector<gap_junction_connection>> connections;
    for (const auto& group: pop) {
        for (auto gid = group.gid_begin; gid < group.gid_end; ++gid) {
            if (!connections.count(gid)) connections.emplace(gid, gn.generate(gid));
        }
    }

    for (const auto& [gid, local_connections]: connections) {
        for (const auto& c: local_connections) {
            EXPECT_TRUE(is_gj_population_member(pop, c.peer));
            EXPECT_TRUE(is_gj_population_member(pop, {gid, c.local}));
            EXPECT_TRUE(selector(c.peer, {gid, c.local}));
            EXPECT_EQ(weight(c.peer, {gid, c.local}), c.weight);
            if (connections.count(c.peer.gid) > 0) {
                const auto& peer_connections = connections[c.peer.gid];
                const auto& g = gid;  // bring into local scope
                if (std::find_if(
                        peer_connections.begin(), peer_connections.end(), [&g, &c](const auto& a) {
                            if (a.peer.gid == g && a.peer.label == c.local) return true;
                            return false;
                        }) == peer_connections.end()) {}
                // make sure connection is generated for peer as well
                EXPECT_TRUE(
                    std::find_if(
                        peer_connections.begin(), peer_connections.end(), [&g, &c](const auto& a) {
                            if (a.peer.gid == g && a.peer.label == c.local) return true;
                            return false;
                        }) != peer_connections.end());
            }
            else {
                const bool symmetric_connections = false;
                EXPECT_TRUE(symmetric_connections);
            }
        }
    }

    return std::accumulate(
        connections.begin(), connections.end(), 0, [](const auto& a, const auto& b) {
            return a + b.second.size();
        });
}

static void for_each_pop_connection(const std::vector<network_cell_group>& pop,
    const std::function<void(const cell_global_label_type&, const cell_global_label_type&)>& func) {

    for (const auto& src: pop) {
        for (const auto& dest: pop) {
            for (auto src_gid = src.gid_begin; src_gid < src.gid_end; ++src_gid) {
                for (auto dest_gid = dest.gid_begin; dest_gid < dest.gid_end; ++dest_gid) {
                    for (const auto& src_label: src.src_labels) {
                        for (const auto& dest_label: dest.dest_labels) {
                            func({src_gid, src_label}, {dest_gid, dest_label});
                        }
                    }
                }
            }
        }
    }
}

static void for_each_pop_connection(const std::vector<spatial_network_cell_group>& pop,
    const std::function<void(const cell_global_label_type&,
        const network_location&,
        const cell_global_label_type&,
        const network_location&)>& func) {

    for (const auto& src: pop) {
        for (const auto& dest: pop) {
            for (auto src_gid = src.gid_begin; src_gid < src.gid_end; ++src_gid) {
                for (auto dest_gid = dest.gid_begin; dest_gid < dest.gid_end; ++dest_gid) {
                    const auto src_location = src.location(src_gid);
                    const auto dest_location = dest.location(dest_gid);

                    for (const auto& src_label: src.src_labels) {
                        for (const auto& dest_label: dest.dest_labels) {
                            func({src_gid, src_label},
                                src_location,
                                {dest_gid, dest_label},
                                dest_location);
                        }
                    }
                }
            }
        }
    }
}

TEST(cell_connection_network, simple) {
    const std::vector<network_cell_group> pop = {{0, 5, {{"a"}}, {}}, {5, 10, {}, {{"b"}}}};
    auto size = check_cell_connections(network_value::uniform(2.0),
        network_value::uniform(3.0),
        network_selection::inter_cell(),
        pop);
    EXPECT_EQ(size, 25);
}

TEST(cell_connection_network, complex) {
    const std::vector<network_cell_group> pop = {
        {0, 5, {{"sa"}, {"sb"}, {"sc"}}, {{"da"}}}, {10, 20, {{"sa"}}, {}}, {20, 50, {}, {{"da"}}}};
    auto size = check_cell_connections(
        network_value(2.0), network_value(3.0), network_selection::all(), pop);
    EXPECT_EQ(size, (3 * 5 + 10) * (5 + 30));
}

TEST(cell_connection_network, no_targets) {
    const std::vector<network_cell_group> pop = {{0, 5, {{"a"}}, {}}, {5, 10, {{"b"}}, {}}};
    auto size = check_cell_connections(
        network_value(2.0), network_value(3.0), network_selection::all(), pop);
    EXPECT_EQ(size, 0);
}

TEST(cell_connection_network, no_sources) {
    const std::vector<network_cell_group> pop = {{0, 5, {}, {{"a"}}}, {5, 10, {}, {{"b"}}}};
    auto size = check_cell_connections(
        network_value(2.0), network_value(3.0), network_selection::all(), pop);
    EXPECT_EQ(size, 0);
}

TEST(cell_connection_network, duplicates) {
    const std::vector<network_cell_group> pop = {{0, 5, {}, {{"a"}}}, {2, 10, {}, {{"b"}}}};
    EXPECT_ANY_THROW(check_cell_connections(
        network_value(2.0), network_value(3.0), network_selection::all(), pop));
}

TEST(gap_junction_network, simple) {
    const std::vector<network_gj_group> pop = {{0, 5, {{"a"}}}, {5, 10, {{"b"}}}};
    auto size = check_gap_junctions(network_value(2.0), network_selection::inter_cell(), pop);
    EXPECT_EQ(size, 100 - 10);
}

TEST(gap_junction_network, complex) {
    const std::vector<network_gj_group> pop = {
        {0, 4, {{"a"}, {"b"}, {"c"}}}, {4, 8, {}}, {10, 25, {{"d"}}}};
    auto size = check_gap_junctions(network_value(2.0), network_selection::inter_cell(), pop);
    EXPECT_EQ(size, (3 * 4 + 15) * (3 * 4 + 15) - 3 * 3 * 4 - 15);  // Substract connections to self
}

TEST(gap_junction_network, empty) {
    const std::vector<network_gj_group> pop = {};
    auto size = check_gap_junctions(network_value(2.0), network_selection::inter_cell(), pop);
    EXPECT_EQ(size, 0);
}

TEST(gap_junction_network, duplicates) {
    const std::vector<network_gj_group> pop = {{0, 5, {{"a"}}}, {2, 10, {{"b"}}}};
    EXPECT_ANY_THROW(check_gap_junctions(network_value(2.0), network_selection::all(), pop));
}

TEST(network_selection, bernoulli_random_100_prob) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::bernoulli_random(42, 1.0);

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_TRUE(s(src, dest));
        });
}

TEST(network_selection, bernoulli_random_0_prob) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::bernoulli_random(42, 0.0);

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_FALSE(s(src, dest));
        });
}

TEST(network_selection, bernoulli_random_consistency) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::bernoulli_random(42, 0.5);

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_EQ(s(src, dest), s(src, dest));
            EXPECT_EQ(s(src, dest), s(dest, src));
        });
}

TEST(network_selection, bernoulli_random_seed) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s1 = network_selection::bernoulli_random(10, 0.5);
    auto s2 = network_selection::bernoulli_random(20, 0.5);

    bool equal_results = true;

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            if (s1(src, dest) != s2(src, dest)) equal_results = false;
        });

    EXPECT_FALSE(equal_results);
}

TEST(network_selection, bernoulli_random_reproducibility) {
    const std::vector<network_cell_group> pop = {{0, 3, {{"a"}}, {}},
        {3, 5, {{"a"}}, {{"d"}}},
        {5, 6, {}, {{"d"}}},
        {7, 11, {{"b"}}, {}},
        {18, 20, {}, {{"e"}}},
        {20, 22, {{"c"}}, {{"e"}}},
        {22, 24, {{"c"}}, {}}};

    auto s = network_selection::bernoulli_random(42, 0.5);
    std::map<cell_global_label_type, std::map<cell_global_label_type, bool>> results;
    results[{0, "a"}][{3, "d"}] = 0;
    results[{0, "a"}][{4, "d"}] = 0;
    results[{0, "a"}][{5, "d"}] = 1;
    results[{1, "a"}][{3, "d"}] = 0;
    results[{1, "a"}][{4, "d"}] = 1;
    results[{1, "a"}][{5, "d"}] = 0;
    results[{2, "a"}][{3, "d"}] = 1;
    results[{2, "a"}][{4, "d"}] = 0;
    results[{2, "a"}][{5, "d"}] = 1;
    results[{3, "a"}][{3, "d"}] = 1;
    results[{3, "a"}][{4, "d"}] = 0;
    results[{3, "a"}][{5, "d"}] = 1;
    results[{4, "a"}][{3, "d"}] = 0;
    results[{4, "a"}][{4, "d"}] = 1;
    results[{4, "a"}][{5, "d"}] = 1;
    results[{0, "a"}][{18, "e"}] = 0;
    results[{0, "a"}][{19, "e"}] = 0;
    results[{0, "a"}][{20, "e"}] = 0;
    results[{0, "a"}][{21, "e"}] = 1;
    results[{1, "a"}][{18, "e"}] = 0;
    results[{1, "a"}][{19, "e"}] = 0;
    results[{1, "a"}][{20, "e"}] = 1;
    results[{1, "a"}][{21, "e"}] = 0;
    results[{2, "a"}][{18, "e"}] = 1;
    results[{2, "a"}][{19, "e"}] = 1;
    results[{2, "a"}][{20, "e"}] = 1;
    results[{2, "a"}][{21, "e"}] = 0;
    results[{3, "a"}][{18, "e"}] = 0;
    results[{3, "a"}][{19, "e"}] = 0;
    results[{3, "a"}][{20, "e"}] = 0;
    results[{3, "a"}][{21, "e"}] = 1;
    results[{4, "a"}][{18, "e"}] = 0;
    results[{4, "a"}][{19, "e"}] = 1;
    results[{4, "a"}][{20, "e"}] = 1;
    results[{4, "a"}][{21, "e"}] = 0;
    results[{7, "b"}][{3, "d"}] = 1;
    results[{7, "b"}][{4, "d"}] = 1;
    results[{7, "b"}][{5, "d"}] = 1;
    results[{8, "b"}][{3, "d"}] = 1;
    results[{8, "b"}][{4, "d"}] = 0;
    results[{8, "b"}][{5, "d"}] = 1;
    results[{9, "b"}][{3, "d"}] = 0;
    results[{9, "b"}][{4, "d"}] = 0;
    results[{9, "b"}][{5, "d"}] = 1;
    results[{10, "b"}][{3, "d"}] = 0;
    results[{10, "b"}][{4, "d"}] = 0;
    results[{10, "b"}][{5, "d"}] = 0;
    results[{7, "b"}][{18, "e"}] = 1;
    results[{7, "b"}][{19, "e"}] = 1;
    results[{7, "b"}][{20, "e"}] = 0;
    results[{7, "b"}][{21, "e"}] = 0;
    results[{8, "b"}][{18, "e"}] = 1;
    results[{8, "b"}][{19, "e"}] = 1;
    results[{8, "b"}][{20, "e"}] = 1;
    results[{8, "b"}][{21, "e"}] = 1;
    results[{9, "b"}][{18, "e"}] = 0;
    results[{9, "b"}][{19, "e"}] = 1;
    results[{9, "b"}][{20, "e"}] = 0;
    results[{9, "b"}][{21, "e"}] = 0;
    results[{10, "b"}][{18, "e"}] = 0;
    results[{10, "b"}][{19, "e"}] = 0;
    results[{10, "b"}][{20, "e"}] = 1;
    results[{10, "b"}][{21, "e"}] = 0;
    results[{20, "c"}][{3, "d"}] = 0;
    results[{20, "c"}][{4, "d"}] = 0;
    results[{20, "c"}][{5, "d"}] = 0;
    results[{21, "c"}][{3, "d"}] = 1;
    results[{21, "c"}][{4, "d"}] = 1;
    results[{21, "c"}][{5, "d"}] = 1;
    results[{22, "c"}][{3, "d"}] = 1;
    results[{22, "c"}][{4, "d"}] = 0;
    results[{22, "c"}][{5, "d"}] = 1;
    results[{23, "c"}][{3, "d"}] = 0;
    results[{23, "c"}][{4, "d"}] = 0;
    results[{23, "c"}][{5, "d"}] = 0;
    results[{20, "c"}][{18, "e"}] = 0;
    results[{20, "c"}][{19, "e"}] = 1;
    results[{20, "c"}][{20, "e"}] = 1;
    results[{20, "c"}][{21, "e"}] = 1;
    results[{21, "c"}][{18, "e"}] = 1;
    results[{21, "c"}][{19, "e"}] = 0;
    results[{21, "c"}][{20, "e"}] = 1;
    results[{21, "c"}][{21, "e"}] = 1;
    results[{22, "c"}][{18, "e"}] = 0;
    results[{22, "c"}][{19, "e"}] = 1;
    results[{22, "c"}][{20, "e"}] = 0;
    results[{22, "c"}][{21, "e"}] = 1;
    results[{23, "c"}][{18, "e"}] = 1;
    results[{23, "c"}][{19, "e"}] = 0;
    results[{23, "c"}][{20, "e"}] = 0;
    results[{23, "c"}][{21, "e"}] = 1;

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            ASSERT_TRUE(results.count(src));
            ASSERT_TRUE(results[src].count(dest));
            EXPECT_EQ(results[src][dest], s(src, dest));
        });
}

TEST(network_selection, custom) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            auto s = network_selection::custom(
                [&](const cell_global_label_type& src_arg, const cell_global_label_type& dest_arg) {
                    EXPECT_EQ(src.gid, src_arg.gid);
                    EXPECT_EQ(dest.gid, dest_arg.gid);
                    EXPECT_EQ(src.label, src_arg.label);
                    EXPECT_EQ(dest.label, dest_arg.label);

                    return true;
                });
            EXPECT_TRUE(s({src.gid, src.label}, {dest.gid, dest.label}));
        });
}

TEST(network_selection, all) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::all();

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_TRUE(s(src, dest));
        });
}

TEST(network_selection, none) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::none();

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_FALSE(s(src, dest));
        });
}

TEST(network_selection, invert) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::invert(network_selection::all());

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_FALSE(s(src, dest));
        });
}

TEST(network_selection, inter_cell) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::inter_cell();

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_EQ(src.gid != dest.gid, s(src, dest));
        });
}

TEST(network_selection, not_equal) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto s = network_selection::not_equal();

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_EQ(src != dest, s(src, dest));
        });
}

TEST(network_selection, and) {
    const cell_global_label_type gl(0, "a");

    const auto true_selection = network_selection::all();
    const auto false_selection = network_selection::invert(network_selection::all());

    EXPECT_TRUE((true_selection & true_selection)(gl, gl));
    EXPECT_FALSE((true_selection & false_selection)(gl, gl));
    EXPECT_FALSE((false_selection & true_selection)(gl, gl));
    EXPECT_FALSE((false_selection & false_selection)(gl, gl));
}

TEST(network_selection, or) {
    const cell_global_label_type gl(0, "a");

    const auto true_selection = network_selection::all();
    const auto false_selection = network_selection::invert(network_selection::all());

    EXPECT_TRUE((true_selection | true_selection)(gl, gl));
    EXPECT_TRUE((true_selection | false_selection)(gl, gl));
    EXPECT_TRUE((false_selection | true_selection)(gl, gl));
    EXPECT_FALSE((false_selection | false_selection)(gl, gl));
}

TEST(network_selection, xor) {
    const cell_global_label_type gl(0, "a");

    const auto true_selection = network_selection::all();
    const auto false_selection = network_selection::invert(network_selection::all());

    EXPECT_FALSE((true_selection ^ true_selection)(gl, gl));
    EXPECT_TRUE((true_selection ^ false_selection)(gl, gl));
    EXPECT_TRUE((false_selection ^ true_selection)(gl, gl));
    EXPECT_FALSE((false_selection ^ false_selection)(gl, gl));
}

TEST(network_value, uniform_distribution) {
    double mean = 0.0;
    std::size_t count = 0;

    const std::vector<network_cell_group> pop = {
        {0, 200, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {250, 270, {{"sc"}}, {{"dc"}}}};

    auto v = network_value::uniform_distribution(42, {-5.0, 3.0});

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            const auto result = v(src, dest);
            mean += result;
            ++count;
            EXPECT_LE(result, 3.0);
            EXPECT_GT(result, -5.0);
        });

    mean /= count;
    EXPECT_NEAR(mean, -1.0, 1e-2);
}

TEST(network_value, uniform_distribution_consistency) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    auto v = network_value::uniform_distribution(42, {2.0, 5.0});

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_EQ(v(src, dest), v(src, dest));
            EXPECT_EQ(v(src, dest), v(dest, src));
        });
}

TEST(network_value, uniform_distribution_reproducibility) {
    const std::vector<network_cell_group> pop = {{0, 2, {{"a"}}, {}},
        {2, 4, {{"a"}}, {{"d"}}},
        {4, 5, {}, {{"d"}}},
        {12, 13, {{"b"}}, {}},
        {13, 14, {{"b"}}, {{"e"}}},
        {14, 15, {}, {{"e"}}}};

    std::map<cell_global_label_type, std::map<cell_global_label_type, double>> results;
    results[{0, "a"}][{2, "d"}] = 0.41817336;
    results[{0, "a"}][{3, "d"}] = 0.27663240;
    results[{0, "a"}][{4, "d"}] = 0.61471849;
    results[{1, "a"}][{2, "d"}] = 0.18429567;
    results[{1, "a"}][{3, "d"}] = 0.19004659;
    results[{1, "a"}][{4, "d"}] = 0.42440975;
    results[{2, "a"}][{2, "d"}] = 0.64354361;
    results[{2, "a"}][{3, "d"}] = 0.82749275;
    results[{2, "a"}][{4, "d"}] = 0.09748023;
    results[{3, "a"}][{2, "d"}] = 0.82380613;
    results[{3, "a"}][{3, "d"}] = 0.64042620;
    results[{3, "a"}][{4, "d"}] = 0.27909534;
    results[{0, "a"}][{13, "e"}] = 0.24041800;
    results[{0, "a"}][{14, "e"}] = 0.55571337;
    results[{1, "a"}][{13, "e"}] = 0.98281076;
    results[{1, "a"}][{14, "e"}] = 0.29195228;
    results[{2, "a"}][{13, "e"}] = 0.84132640;
    results[{2, "a"}][{14, "e"}] = 0.21757594;
    results[{3, "a"}][{13, "e"}] = 0.48160950;
    results[{3, "a"}][{14, "e"}] = 0.12191449;
    results[{12, "b"}][{2, "d"}] = 0.41704783;
    results[{12, "b"}][{3, "d"}] = 0.84197325;
    results[{12, "b"}][{4, "d"}] = 0.47905073;
    results[{13, "b"}][{2, "d"}] = 0.19514098;
    results[{13, "b"}][{3, "d"}] = 0.01418517;
    results[{13, "b"}][{4, "d"}] = 0.47480690;
    results[{12, "b"}][{13, "e"}] = 0.21116588;
    results[{12, "b"}][{14, "e"}] = 0.79147732;
    results[{13, "b"}][{13, "e"}] = 0.48435615;
    results[{13, "b"}][{14, "e"}] = 0.55361571;

    auto v = network_value::uniform_distribution(42, {0.0, 1.0});

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            ASSERT_TRUE(results.count(src));
            ASSERT_TRUE(results[src].count(dest));
            EXPECT_NEAR(v(src, dest), results[src][dest], 1e-7);
        });
}

TEST(network_value, normal_distribution) {
    const std::vector<network_cell_group> pop = {{0, 500, {{"a"}, {"b"}}, {{"d"}}}};

    const double mean = 5.0;
    const double std_dev = 3.0;

    auto v = network_value::normal_distribution(42, mean, std_dev);

    double sample_mean = 0.0;
    double sample_dev = 0.0;

    std::size_t count = 0;

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            const auto result = v(src, dest);
            sample_mean += result;
            sample_dev += (result - mean) * (result - mean);
            ++count;
        });

    sample_mean /= count;
    sample_dev = std::sqrt(sample_dev / count);

    EXPECT_NEAR(sample_mean, mean, 1e-2);
    EXPECT_NEAR(sample_dev, std_dev, 1e-2);
}

TEST(network_value, normal_distribution_reproducibility) {
    const std::vector<network_cell_group> pop = {{0, 2, {{"a"}}, {}},
        {2, 4, {{"a"}}, {{"d"}}},
        {4, 5, {}, {{"d"}}},
        {12, 13, {{"b"}}, {}},
        {13, 14, {{"b"}}, {{"e"}}},
        {14, 15, {}, {{"e"}}}};

    std::map<cell_global_label_type, std::map<cell_global_label_type, double>> results;
    results[{0, "a"}][{2, "d"}] = 0.08772826;
    results[{0, "a"}][{3, "d"}] = 9.67308342;
    results[{0, "a"}][{4, "d"}] = 3.73949064;
    results[{1, "a"}][{2, "d"}] = -0.56031250;
    results[{1, "a"}][{3, "d"}] = 2.37198282;
    results[{1, "a"}][{4, "d"}] = 6.14574844;
    results[{2, "a"}][{2, "d"}] = -1.51583318;
    results[{2, "a"}][{3, "d"}] = 7.05109701;
    results[{2, "a"}][{4, "d"}] = 2.35098292;
    results[{3, "a"}][{2, "d"}] = 5.55626489;
    results[{3, "a"}][{3, "d"}] = 8.07950315;
    results[{3, "a"}][{4, "d"}] = 5.53472845;
    results[{0, "a"}][{13, "e"}] = 4.72865713;
    results[{0, "a"}][{14, "e"}] = 5.91688185;
    results[{1, "a"}][{13, "e"}] = 10.81253890;
    results[{1, "a"}][{14, "e"}] = 8.86701434;
    results[{2, "a"}][{13, "e"}] = 1.16593915;
    results[{2, "a"}][{14, "e"}] = 5.09421785;
    results[{3, "a"}][{13, "e"}] = 4.76063628;
    results[{3, "a"}][{14, "e"}] = 1.69353474;
    results[{12, "b"}][{2, "d"}] = 4.76564698;
    results[{12, "b"}][{3, "d"}] = 5.81550403;
    results[{12, "b"}][{4, "d"}] = 13.99375691;
    results[{13, "b"}][{2, "d"}] = 6.03529503;
    results[{13, "b"}][{3, "d"}] = 1.96907188;
    results[{13, "b"}][{4, "d"}] = 5.38927382;
    results[{12, "b"}][{13, "e"}] = 5.89951548;
    results[{12, "b"}][{14, "e"}] = 3.33922968;
    results[{13, "b"}][{13, "e"}] = 12.49977574;
    results[{13, "b"}][{14, "e"}] = 3.97665339;

    const double mean = 5.0;
    const double std_dev = 3.0;

    auto v = network_value::normal_distribution(42, mean, std_dev);

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            ASSERT_TRUE(results.count(src));
            ASSERT_TRUE(results[src].count(dest));
            EXPECT_NEAR(v(src, dest), results[src][dest], 1e-7);
        });
}

TEST(network_value, truncated_normal_distribution) {
    const std::vector<network_cell_group> pop = {{0, 500, {{"a"}, {"b"}}, {{"d"}}}};

    const double mean = 5.0;
    const double std_dev = 3.0;
    const double lower_bound = 1.0;
    const double upper_bound = 9.0;

    auto v =
        network_value::truncated_normal_distribution(42, mean, std_dev, {lower_bound, upper_bound});

    double sample_mean = 0.0;

    std::size_t count = 0;

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            const auto result = v(src, dest);
            EXPECT_GT(result, lower_bound);
            EXPECT_LE(result, upper_bound);
            sample_mean += result;
            ++count;
        });

    sample_mean /= count;

    EXPECT_NEAR(sample_mean, mean, 1e-2);
}

TEST(network_value, truncated_normal_distribution_reproducibility) {
    const std::vector<network_cell_group> pop = {{0, 2, {{"a"}}, {}},
        {2, 4, {{"a"}}, {{"d"}}},
        {4, 5, {}, {{"d"}}},
        {12, 13, {{"b"}}, {}},
        {13, 14, {{"b"}}, {{"e"}}},
        {14, 15, {}, {{"e"}}}};

    std::map<cell_global_label_type, std::map<cell_global_label_type, double>> results;
    results[{0, "a"}][{2, "d"}] = 6.75583032;
    results[{0, "a"}][{3, "d"}] = 4.45202654;
    results[{0, "a"}][{4, "d"}] = 7.94629345;
    results[{1, "a"}][{2, "d"}] = 5.38361275;
    results[{1, "a"}][{3, "d"}] = 4.66503239;
    results[{1, "a"}][{4, "d"}] = 5.52403259;
    results[{2, "a"}][{2, "d"}] = 3.61250351;
    results[{2, "a"}][{3, "d"}] = 6.23056041;
    results[{2, "a"}][{4, "d"}] = 5.90644730;
    results[{3, "a"}][{2, "d"}] = 6.14728960;
    results[{3, "a"}][{3, "d"}] = 5.22867021;
    results[{3, "a"}][{4, "d"}] = 4.42161676;
    results[{0, "a"}][{13, "e"}] = 5.99341261;
    results[{0, "a"}][{14, "e"}] = 3.37764114;
    results[{1, "a"}][{13, "e"}] = 7.27707923;
    results[{1, "a"}][{14, "e"}] = 4.00826224;
    results[{2, "a"}][{13, "e"}] = 6.64256804;
    results[{2, "a"}][{14, "e"}] = 5.20821494;
    results[{3, "a"}][{13, "e"}] = 7.06933652;
    results[{3, "a"}][{14, "e"}] = 5.91076864;
    results[{12, "b"}][{2, "d"}] = 4.70123400;
    results[{12, "b"}][{3, "d"}] = 4.80106583;
    results[{12, "b"}][{4, "d"}] = 3.47787931;
    results[{13, "b"}][{2, "d"}] = 4.64732400;
    results[{13, "b"}][{3, "d"}] = 5.81092091;
    results[{13, "b"}][{4, "d"}] = 3.49918183;
    results[{12, "b"}][{13, "e"}] = 7.75582442;
    results[{12, "b"}][{14, "e"}] = 4.20438519;
    results[{13, "b"}][{13, "e"}] = 5.40959364;
    results[{13, "b"}][{14, "e"}] = 5.16618053;

    const double mean = 5.0;
    const double std_dev = 3.0;

    auto v =
        network_value::truncated_normal_distribution(42, mean, std_dev, {mean - 2.0, mean + 3.0});

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            ASSERT_TRUE(results.count(src));
            ASSERT_TRUE(results[src].count(dest));
            EXPECT_NEAR(v(src, dest), results[src][dest], 1e-7);
        });
}

TEST(network_value, custom) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            auto s = network_value::custom(
                [&](const cell_global_label_type& src_arg, const cell_global_label_type& dest_arg) {
                    EXPECT_EQ(src.gid, src_arg.gid);
                    EXPECT_EQ(dest.gid, dest_arg.gid);
                    EXPECT_EQ(src.label, src_arg.label);
                    EXPECT_EQ(dest.label, dest_arg.label);

                    return 2.0;
                });

            EXPECT_DOUBLE_EQ(s(src, dest), 2.0);
        });
}

TEST(network_value, uniform) {
    const std::vector<network_cell_group> pop = {
        {0, 20, {{"sa"}, {"sb"}}, {{"da"}, {"db"}}}, {25, 42, {{"sc"}}, {}}};

    const auto v1 = network_value::uniform(5.0);
    const network_value v2(5.0);

    for_each_pop_connection(
        pop, [&](const cell_global_label_type& src, const cell_global_label_type& dest) {
            EXPECT_DOUBLE_EQ(v1(src, dest), 5.0);
            EXPECT_DOUBLE_EQ(v1(src, dest), v2(src, dest));
        });
}

TEST(spatial_network_value, conversion) {
    const network_location loc = {0.0, 0.0, 0.0};
    const std::vector<spatial_network_cell_group> pop = {
        {0, {{"a"}}, {{"b"}}, std::vector<network_location>(50, loc)}};

    const auto v1 = network_value::normal_distribution(42, 2.0, 2.0);
    const spatial_network_value v2(v1);

    for_each_pop_connection(pop,
        [&](const cell_global_label_type& src,
            const network_location& src_location,
            const cell_global_label_type& dest,
            const network_location& dest_location) {
            EXPECT_DOUBLE_EQ(v1(src, dest), v2(src, src_location, dest, dest_location));
        });
}


TEST(spatial_network_value, custom) {

    const std::vector<spatial_network_cell_group> pop = {
        {0, {{"a"}}, {{"b"}}, std::vector<network_location>(10, network_location{1.0, 1.0, 1.0})},
        {10, {{"c"}}, {{"d"}}, std::vector<network_location>(10, network_location{2.0, 2.0, 2.0})}};

    for_each_pop_connection(pop,
        [&](const cell_global_label_type& src,
            const network_location& src_location,
            const cell_global_label_type& dest,
            const network_location& dest_location) {
            auto s = spatial_network_value::custom([&](const cell_global_label_type& src_arg,
                                                       const network_location& src_location_arg,
                                                       const cell_global_label_type& dest_arg,
                                                       const network_location& dest_location_arg,
                                                       double distance_arg) {
                EXPECT_EQ(src.gid, src_arg.gid);
                EXPECT_EQ(dest.gid, dest_arg.gid);
                EXPECT_EQ(src.label, src_arg.label);
                EXPECT_EQ(dest.label, dest_arg.label);
                EXPECT_EQ(dest_location, src_location);

                return 2.0;
            });

            EXPECT_DOUBLE_EQ(s(src, src_location, dest, dest_location), 2.0);
        });
}
