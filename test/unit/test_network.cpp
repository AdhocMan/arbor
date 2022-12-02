#include <gtest/gtest.h>

#include <arbor/network.hpp>

#include <algorithm>
#include <numeric>
#include <unordered_map>

using namespace arb;

static bool is_population_member(const network_population& pop, const cell_global_label_type& l) {
    for (const auto& range: pop) {
        if (range.begin <= l.gid && range.end > l.gid && range.label == l.label) return true;
    }
    return false;
}

static long check_cell_connections(network_value weight,
    network_value delay,
    network_selection selector,
    const network_population& src_pop,
    const network_population& dest_pop) {

    cell_connection_network cn(weight, delay, selector, src_pop, dest_pop);

    // A gid may occur more than once, so gather unique sizes first
    std::unordered_map<cell_gid_type, std::size_t> sizes;

    for (const auto& dest: dest_pop) {
        for (auto gid = dest.begin; gid < dest.end; ++gid) {
            auto connections = cn.generate(gid);
            sizes.insert_or_assign(gid, connections.size());
            for (const auto& c: connections) {
                EXPECT_TRUE(is_population_member(src_pop, c.source));
                EXPECT_TRUE(is_population_member(dest_pop, {gid, c.dest}));
                EXPECT_TRUE(selector(c.source, {gid, c.dest}));
                EXPECT_EQ(weight(c.source, {gid, c.dest}), c.weight);
                EXPECT_EQ(delay(c.source, {gid, c.dest}), c.delay);
            }
        }
    }

    return std::accumulate(
        sizes.begin(), sizes.end(), 0, [](const auto& a, const auto& b) { return a + b.second; });
}

static long check_gap_junctions(network_value weight,
    network_selection selector,
    const network_population& src_pop,
    const network_population& dest_pop) {

    gap_junction_network gn(weight, selector, src_pop, dest_pop);

    std::unordered_map<cell_gid_type, std::vector<gap_junction_connection>> connections;
    for (const auto& dest: dest_pop) {
        for (auto gid = dest.begin; gid < dest.end; ++gid) {
            if (!connections.count(gid)) connections.emplace(gid, gn.generate(gid));
        }
    }
    for (const auto& src: src_pop) {
        for (auto gid = src.begin; gid < src.end; ++gid) {
            if (!connections.count(gid)) connections.emplace(gid, gn.generate(gid));
        }
    }

    for (const auto& [gid, local_connections]: connections) {
        for (const auto& c: local_connections) {
            EXPECT_TRUE(
                is_population_member(src_pop, c.peer) || is_population_member(dest_pop, c.peer));
            EXPECT_TRUE(is_population_member(dest_pop, {gid, c.local}) ||
                        is_population_member(src_pop, {gid, c.local}));
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

TEST(network_population, unique) {
    const network_population pop = {{10, 15, "a"}, {5, 12, "a"}, {13, 20, "a"}, {13, 20, "b"}};
    auto unique_pop = unique(pop);
    network_population reference_pop = {{5, 20, "a"}, {13, 20, "b"}};
    std::sort(unique_pop.begin(), unique_pop.end());
    std::sort(reference_pop.begin(), reference_pop.end());

    ASSERT_EQ(unique_pop.size(), reference_pop.size());
    for (std::size_t i = 0; i < reference_pop.size(); ++i) {
        EXPECT_EQ(unique_pop[i], reference_pop[i]);
    }
}

TEST(network_population, join_two) {
    const network_population pop1 = {{10, 15, "a"}, {5, 12, "a"}};
    const network_population pop2 = {{13, 20, "a"}, {13, 20, "b"}};
    network_population reference_pop = {{10, 15, "a"}, {5, 12, "a"}, {13, 20, "a"}, {13, 20, "b"}};
    auto joined_pop = join(pop1, pop2);
    std::sort(joined_pop.begin(), joined_pop.end());
    std::sort(reference_pop.begin(), reference_pop.end());

    ASSERT_EQ(joined_pop.size(), reference_pop.size());
    for (std::size_t i = 0; i < reference_pop.size(); ++i) {
        EXPECT_EQ(joined_pop[i], reference_pop[i]);
    }
}

TEST(network_population, join_multiple) {
    const network_population pop1 = {{10, 15, "a"}, {5, 12, "a"}};
    const network_population pop2 = {{13, 20, "a"}, {13, 20, "b"}};
    const network_population pop3 = {{0, 5, "c"}};
    const network_population pop4 = {{40, 50, "d"}, {60, 65, "b"}};
    network_population reference_pop = {{10, 15, "a"},
        {5, 12, "a"},
        {13, 20, "a"},
        {13, 20, "b"},
        {0, 5, "c"},
        {40, 50, "d"},
        {60, 65, "b"}};
    auto joined_pop = join(pop1, pop2, pop3, pop4);
    std::sort(joined_pop.begin(), joined_pop.end());
    std::sort(reference_pop.begin(), reference_pop.end());

    ASSERT_EQ(joined_pop.size(), reference_pop.size());
    for (std::size_t i = 0; i < reference_pop.size(); ++i) {
        EXPECT_EQ(joined_pop[i], reference_pop[i]);
    }
}

TEST(cell_connection_network, simple) {
    const network_population src_pop = {{0, 5, "a"}};
    const network_population dest_pop = {{5, 10, "b"}};
    auto size =
        check_cell_connections(2.0, 3.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 25);
}

TEST(cell_connection_network, complex) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};
    auto size =
        check_cell_connections(2.0, 3.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, (3 + 37) * (5 + 8 + 130) - (2 + 2 + 30));  // Substract connections to self
}

TEST(cell_connection_network, no_targets) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {};
    auto size =
        check_cell_connections(2.0, 3.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 0);
}

TEST(cell_connection_network, no_sources) {
    const network_population src_pop = {};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};
    auto size =
        check_cell_connections(2.0, 3.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 0);
}

TEST(cell_connection_network, duplicates) {
    const network_population src_pop = {{0, 10, "a"}, {0, 10, "a"}};
    const network_population dest_pop = {{0, 10, "b"}, {0, 10, "b"}};
    auto size = check_cell_connections(2.0, 3.0, network_selection::all(), src_pop, dest_pop);
    EXPECT_EQ(size, 100);
}

TEST(gap_junction_network, simple) {
    const network_population src_pop = {{0, 5, "a"}};
    const network_population dest_pop = {{5, 10, "b"}};
    auto size = check_gap_junctions(2.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 2 * 25);
}

TEST(gap_junction_network, complex) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};
    auto size = check_gap_junctions(2.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(
        size, 2 * ((3 + 37) * (5 + 8 + 130) - (2 + 2 + 30)));  // Substract connections to self
}

TEST(gap_junction_network, empty) {
    const network_population src_pop = {};
    const network_population dest_pop = {};
    auto size = check_gap_junctions(2.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 0);
}

TEST(gap_junction_network, no_sources) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {};
    auto size = check_gap_junctions(2.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 0);
}

TEST(gap_junction_network, no_targets) {
    const network_population src_pop = {};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};
    auto size = check_gap_junctions(2.0, network_selection::inter_cell(), src_pop, dest_pop);
    EXPECT_EQ(size, 0);
}

TEST(gap_junction_network, duplicates) {
    const network_population src_pop = {{0, 10, "a"}};
    const network_population dest_pop = {{0, 10, "a"}};
    auto size = check_gap_junctions(2.0, network_selection::all(), src_pop, dest_pop);
    EXPECT_EQ(size, 100);
}

TEST(network_selection, bernoulli_random_100_prob) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s = network_selection::bernoulli_random(42, 1.0);
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    EXPECT_TRUE(s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, bernoulli_random_0_prob) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s = network_selection::bernoulli_random(42, 0.0);
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    EXPECT_FALSE(s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, bernoulli_random_consistency) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s = network_selection::bernoulli_random(42, 0.5);
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);
                    EXPECT_EQ(s(src_gl, dest_gl), s(src_gl, dest_gl));
                    EXPECT_EQ(s(src_gl, dest_gl), s(dest_gl, src_gl));
                }
            }
        }
    }
}

TEST(network_selection, bernoulli_random_seed) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s1 = network_selection::bernoulli_random(10, 0.5);
    auto s2 = network_selection::bernoulli_random(20, 0.5);

    bool equal_results = true;
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);
                    if (s1(src_gl, dest_gl) != s2(src_gl, dest_gl)) equal_results = false;
                }
            }
        }
    }

    EXPECT_FALSE(equal_results);
}

TEST(network_selection, bernoulli_random_reproducibility) {
    const network_population src_pop = {{0, 5, "a"}, {7, 11, "b"}, {20, 24, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {18, 22, "e"}};

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

    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);
                    ASSERT_TRUE(results.count(src_gl));
                    ASSERT_TRUE(results[src_gl].count(dest_gl));
                    EXPECT_EQ(results[src_gl][dest_gl], s(src_gl, dest_gl));
                }
            }
        }
    }
}

TEST(network_selection, custom) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {

                    auto s = network_selection::custom([&](const cell_global_label_type& src_arg,
                                                           const cell_global_label_type& dest_arg) {
                        EXPECT_EQ(src_gid, src_arg.gid);
                        EXPECT_EQ(dest_gid, dest_arg.gid);
                        EXPECT_EQ(src.label, src_arg.label);
                        EXPECT_EQ(dest.label, dest_arg.label);

                        return true;
                    });
                    EXPECT_TRUE(s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, all) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s = network_selection::all();
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    EXPECT_TRUE(s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, none) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s = network_selection::none();
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    EXPECT_FALSE(s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, invert) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto s = network_selection::invert(network_selection::all());
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    EXPECT_FALSE(s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, inter_cell) {
    const network_population src_pop = {{0, 5, "a"}, {10, 25, "b"}};
    const network_population dest_pop = {{3, 8, "a"}, {20, 30, "e"}};

    auto s = network_selection::inter_cell();
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    EXPECT_EQ(src_gid != dest_gid, s({src_gid, src.label}, {dest_gid, dest.label}));
                }
            }
        }
    }
}

TEST(network_selection, not_equal) {
    const network_population src_pop = {{0, 5, "a"}, {10, 25, "b"}};
    const network_population dest_pop = {{2, 8, "a"}, {20, 30, "e"}};

    auto s = network_selection::not_equal();
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);

                    EXPECT_EQ(src_gl != dest_gl, s(src_gl, dest_gl));
                }
            }
        }
    }
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
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto v = network_value::uniform_distribution(42, {-5.0, 3.0});
    double mean = 0.0;
    std::size_t count = 0;
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid, ++count) {
                    const auto result = v({src_gid, src.label}, {dest_gid, dest.label});
                    mean += result;
                    EXPECT_LE(result, 3.0);
                    EXPECT_GT(result, -5.0);
                }
            }
        }
    }
    mean /= count;
    EXPECT_NEAR(mean, -1.0, 1e-2);
}

TEST(network_value, uniform_distribution_consistency) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    auto v = network_value::uniform_distribution(42, {2.0, 5.0});
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);
                    EXPECT_EQ(v(src_gl, dest_gl), v(src_gl, dest_gl));
                    EXPECT_EQ(v(src_gl, dest_gl), v(dest_gl, src_gl));
                }
            }
        }
    }
}

TEST(network_value, uniform_distribution_reproducibility) {
    const network_population src_pop = {{0, 4, "a"}, {12, 14, "b"}};
    const network_population dest_pop = {{2, 5, "d"}, {13, 15, "e"}};

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
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);
                    ASSERT_TRUE(results.count(src_gl));
                    ASSERT_TRUE(results[src_gl].count(dest_gl));
                    EXPECT_NEAR(v(src_gl, dest_gl), results[src_gl][dest_gl], 1e-7);
                }
            }
        }
    }
}

TEST(network_value, normal_distribution) {
    const network_population src_pop = {{0, 100, "a"}, {100, 500, "b"}};
    const network_population dest_pop = {{0, 500, "d"}};

    const double mean = 5.0;
    const double std_dev = 3.0;

    auto v = network_value::normal_distribution(42, mean, std_dev);

    double sample_mean = 0.0;
    double sample_dev = 0.0;

    std::size_t count = 0;
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid, ++count) {
                    const auto result = v({src_gid, src.label}, {dest_gid, dest.label});
                    sample_mean += result;
                    sample_dev += (result - mean) * (result - mean);
                }
            }
        }
    }

    sample_mean /= count;
    sample_dev = std::sqrt(sample_dev / count);

    EXPECT_NEAR(sample_mean, mean, 1e-2);
    EXPECT_NEAR(sample_dev, std_dev, 1e-2);
}

TEST(network_value, normal_distribution_reproducibility) {
    const network_population src_pop = {{0, 4, "a"}, {12, 14, "b"}};
    const network_population dest_pop = {{2, 5, "d"}, {13, 15, "e"}};

    std::map<cell_global_label_type, std::map<cell_global_label_type, double>> results;
    results[{0, "a"}][{2, "d"}] = 5.22844040;
    results[{0, "a"}][{3, "d"}] = 12.11075575;
    results[{0, "a"}][{4, "d"}] = 1.79567588;
    results[{1, "a"}][{2, "d"}] = 6.40208047;
    results[{1, "a"}][{3, "d"}] = 8.88639034;
    results[{1, "a"}][{4, "d"}] = 5.66168005;
    results[{2, "a"}][{2, "d"}] = 4.71464577;
    results[{2, "a"}][{3, "d"}] = 3.09002168;
    results[{2, "a"}][{4, "d"}] = 8.42688714;
    results[{3, "a"}][{2, "d"}] = 2.81869263;
    results[{3, "a"}][{3, "d"}] = 0.88072245;
    results[{3, "a"}][{4, "d"}] = 11.86901555;
    results[{0, "a"}][{13, "e"}] = 9.06395119;
    results[{0, "a"}][{14, "e"}] = 4.45229310;
    results[{1, "a"}][{13, "e"}] = 4.61727001;
    results[{1, "a"}][{14, "e"}] = 6.82659744;
    results[{2, "a"}][{13, "e"}] = 1.77223483;
    results[{2, "a"}][{14, "e"}] = 12.92401473;
    results[{3, "a"}][{13, "e"}] = 5.56840067;
    results[{3, "a"}][{14, "e"}] = 6.97205634;
    results[{12, "b"}][{2, "d"}] = 6.15835914;
    results[{12, "b"}][{3, "d"}] = 1.41352179;
    results[{12, "b"}][{4, "d"}] = 5.24003553;
    results[{13, "b"}][{2, "d"}] = 7.19132905;
    results[{13, "b"}][{3, "d"}] = 5.20244979;
    results[{13, "b"}][{4, "d"}] = 5.53694682;
    results[{12, "b"}][{13, "e"}] = 6.58797260;
    results[{12, "b"}][{14, "e"}] = -0.74208974;
    results[{13, "b"}][{13, "e"}] = 5.13288963;
    results[{13, "b"}][{14, "e"}] = 4.84014301;

    const double mean = 5.0;
    const double std_dev = 3.0;

    auto v = network_value::normal_distribution(42, mean, std_dev);
    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);

                    ASSERT_TRUE(results.count(src_gl));
                    ASSERT_TRUE(results[src_gl].count(dest_gl));
                    EXPECT_NEAR(v(src_gl, dest_gl), results[src_gl][dest_gl], 1e-7);
                }
            }
        }
    }
}

TEST(network_value, custom) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {

                    auto s = network_value::custom([&](const cell_global_label_type& src_arg,
                                                       const cell_global_label_type& dest_arg) {
                        EXPECT_EQ(src_gid, src_arg.gid);
                        EXPECT_EQ(dest_gid, dest_arg.gid);
                        EXPECT_EQ(src.label, src_arg.label);
                        EXPECT_EQ(dest.label, dest_arg.label);

                        return 2.0;
                    });

                    EXPECT_DOUBLE_EQ(s({src_gid, src.label}, {dest_gid, dest.label}), 2.0);
                }
            }
        }
    }
}

TEST(network_value, uniform) {
    const network_population src_pop = {{0, 5, "a"}, {7, 15, "b"}, {20, 150, "c"}};
    const network_population dest_pop = {{3, 6, "d"}, {13, 50, "e"}};

    const auto v1 = network_value::uniform(5.0);
    const network_value v2(5.0);

    for (const auto& src: src_pop) {
        for (const auto& dest: dest_pop) {
            for (auto src_gid = src.begin; src_gid < src.end; ++src_gid) {
                for (auto dest_gid = dest.begin; dest_gid < dest.end; ++dest_gid) {
                    const auto src_gl = cell_global_label_type(src_gid, src.label);
                    const auto dest_gl = cell_global_label_type(dest_gid, dest.label);

                    EXPECT_DOUBLE_EQ(v1(src_gl, dest_gl), 5.0);
                    EXPECT_DOUBLE_EQ(v1(src_gl, dest_gl), v2(src_gl, dest_gl));
                }
            }
        }
    }
}
