#include "../gtest.h"
#include "arbor/morph/primitives.hpp"

#include <arbor/iexpr.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/cable_cell_param.hpp>

#include <arborio/label_parse.hpp>

using namespace arb;
using namespace arborio::literals;

TEST(iexpr, distance_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;

    // test distance to single point
    auto dist = thingify(arb::iexpr::distance(scale, arb::mlocation{1, 0.2}), prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), scale * 12.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), scale * 22.0);

    // test distance to multiple points
    dist = thingify(arb::iexpr::distance(scale, arb::mlocation_list({arb::mlocation{1, 1.0}, arb::mlocation{2, 1.0}})), prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), scale * 15.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), scale * 5.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), scale * 10.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), scale * 30.0);
}

TEST(iexpr, distance_region) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 3.0;

    // test distance to single cable region
    auto dist = thingify(arb::iexpr::distance(scale, arb::mcable{0, 0.2, 0.8}), prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), scale * 12.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), scale * 12.0);

    // test distance to multi cable region
    dist = thingify(arb::iexpr::distance(scale, arb::mcable_list{arb::mcable{1, 0.2, 0.8}, arb::mcable{2, 0.6, 1.0}}), prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), scale * 2.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), scale * 22.0);
}

TEST(iexpr, proximal_distance_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;

    // test distance to single point
    auto dist = thingify(arb::iexpr::proximal_distance(scale, arb::mlocation{1, 0.2}), prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multiple points
    dist = thingify(arb::iexpr::proximal_distance(scale,
                        arb::mlocation_list({arb::mlocation{1, 0.2}, arb::mlocation{2, 1.0}})),
        prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), scale * 7);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), scale * 10.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), 0.0);
}

TEST(iexpr, distal_distance_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;

    // test distance to single point
    auto dist = thingify(arb::iexpr::distal_distance(scale, arb::mlocation{1, 0.2}), prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multiple points
    dist = thingify(arb::iexpr::distal_distance(scale,
                        arb::mlocation_list({arb::mlocation{1, 0.2}, arb::mlocation{3, 0.4}})),
        prov);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(dist->eval(prov, {3, 0.0, 1.0}), scale * 2.0);
}

TEST(iexpr, scalar) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::scalar(2.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0);
}

TEST(iexpr, radius) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 20, 2}, 2);
    tree.append(0,     {0, 0, 10, 10}, {0, 0, 30, 5}, 3);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;
    auto e = thingify(arb::iexpr::radius(scale), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), scale * 10.0);
    EXPECT_DOUBLE_EQ(e->eval(prov, {1, 0.0, 1.0}), scale * 1.5);
    EXPECT_DOUBLE_EQ(e->eval(prov, {2, 0.0, 1.0}), scale * 7.5);
}

TEST(iexpr, diameter) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0,     {0, 0, 10, 1}, {0, 0, 20, 2}, 2);
    tree.append(0,     {0, 0, 10, 10}, {0, 0, 30, 5}, 3);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;
    auto e = thingify(arb::iexpr::diameter(scale), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), scale * 20.0);
    EXPECT_DOUBLE_EQ(e->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(e->eval(prov, {2, 0.0, 1.0}), scale * 15.0);
}

TEST(iexpr, add) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::add(arb::iexpr::scalar(2.0), arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 + 3.0 * 10.0);
}

TEST(iexpr, mul) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::mul(arb::iexpr::scalar(2.0), arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 * 3.0 * 10.0);
}
