#include "../gtest.h"
#include "arbor/morph/primitives.hpp"

#include <arbor/cable_cell.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/iexpr.hpp>
#include <arborio/label_parse.hpp>

#include <cmath>

using namespace arb;
using namespace arborio::literals;

TEST(iexpr, distance_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;

    // test distance to single point
    auto ex = thingify(arb::iexpr::distance(scale, arb::mlocation{1, 0.2}), prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 12.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), scale * 22.0);

    // test distance to multiple points
    ex = thingify(arb::iexpr::distance(
                      scale, arb::mlocation_list({arb::mlocation{1, 1.0}, arb::mlocation{2, 1.0}})),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 15.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), scale * 5.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 10.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), scale * 30.0);
}

TEST(iexpr, distance_region) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 3.0;

    // test distance to single cable region
    auto ex = thingify(arb::iexpr::distance(scale, arb::mcable{0, 0.2, 0.8}), prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 12.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), scale * 12.0);

    // test distance to multi cable region
    ex = thingify(arb::iexpr::distance(
                      scale, arb::mcable_list{arb::mcable{1, 0.2, 0.8}, arb::mcable{2, 0.6, 1.0}}),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 2.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), scale * 22.0);
}

TEST(iexpr, proximal_distance_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;

    // test distance to single point
    auto ex = thingify(arb::iexpr::proximal_distance(scale, arb::mlocation{1, 0.2}), prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multiple points
    ex = thingify(arb::iexpr::proximal_distance(
                      scale, arb::mlocation_list({arb::mlocation{1, 0.2}, arb::mlocation{2, 1.0}})),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 7);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 10.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);
}

TEST(iexpr, proximal_distance_region) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 3.0;

    // test distance to single cable region
    auto ex = thingify(arb::iexpr::proximal_distance(scale, arb::mcable{1, 0.2, 0.8}), prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multi cable region
    ex = thingify(
        arb::iexpr::proximal_distance(scale,
            arb::mcable_list{
                arb::mcable{1, 0.2, 0.8}, arb::mcable{2, 0.6, 1.0}, arb::mcable{3, 0.1, 0.3}}),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), scale * 7.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 2.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);
}

TEST(iexpr, distal_distance_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 2.0;

    // test distance to single point
    auto ex = thingify(arb::iexpr::distal_distance(scale, arb::mlocation{1, 0.2}), prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multiple points
    ex = thingify(arb::iexpr::distal_distance(
                      scale, arb::mlocation_list({arb::mlocation{1, 0.2}, arb::mlocation{3, 0.4}})),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), scale * 3.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), scale * 2.0);
}

TEST(iexpr, distal_distance_region) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double scale = 3.0;

    // test distance to single cable region
    auto ex = thingify(arb::iexpr::distal_distance(scale, arb::mcable{1, 0.2, 0.8}), prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multi cable region
    ex = thingify(
        arb::iexpr::distal_distance(scale,
            arb::mcable_list{
                arb::mcable{1, 0.2, 0.8}, arb::mcable{2, 0.3, 0.4}, arb::mcable{3, 0.1, 0.2}}),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), scale * 2.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), scale * 6.0);
}

TEST(iexpr, interpolation_locset) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double prox_value = 2.0;
    const double dist_value = 3.0;

    // test single point
    auto ex = thingify(arb::iexpr::interpolation(
                           prox_value, arb::mlocation{1, 0.2}, dist_value, arb::mlocation{1, 0.2}),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test evaluation on ends of interval
    ex = thingify(arb::iexpr::interpolation(
                      prox_value, arb::mlocation{1, 0.2}, dist_value, arb::mlocation{1, 0.8}),
        prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.2, 0.2}), prox_value);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.8, 0.8}), dist_value);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multiple points
    ex =
        thingify(arb::iexpr::interpolation(prox_value,
                     arb::mlocation{0, 0.3},
                     dist_value,
                     arb::mlocation_list(
                         {arb::mlocation{1, 0.2}, arb::mlocation{2, 0.8}, arb::mlocation{3, 1.0}})),
            prov);
    EXPECT_DOUBLE_EQ(
        ex->eval(prov, {0, 0.0, 1.0}), 7.0 / 9.0 * prox_value + 2.0 / 9.0 * dist_value);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(
        ex->eval(prov, {2, 0.0, 1.0}), 6.0 / 23.0 * prox_value + 17.0 / 23.0 * dist_value);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);
}

TEST(iexpr, interpolation_region) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 1}, 3);
    tree.append(0, {0, 0, 10, 1}, {0, 0, 30, 1}, 4);
    tree.append(mnpos, {0, 0, 0, 2}, {0, 0, -20, 2}, 2);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    const double prox_value = 2.0;
    const double dist_value = 3.0;

    // test distance to single point
    auto ex =
        thingify(arb::iexpr::interpolation(
                     prox_value, arb::mcable{1, 0.2, 0.7}, dist_value, arb::mcable{1, 0.2, 0.7}),
            prov);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {0, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);

    // test distance to multiple points
    ex = thingify(
        arb::iexpr::interpolation(prox_value,
            arb::mcable{0, 0.1, 0.3},
            dist_value,
            arb::mcable_list(
                {arb::mcable{1, 0.2, 0.4}, arb::mcable{2, 0.2, 1.0}, arb::mcable{3, 0.1, 1.0}})),
        prov);
    EXPECT_DOUBLE_EQ(
        ex->eval(prov, {0, 0.0, 1.0}), 7.0 / 9.0 * prox_value + 2.0 / 9.0 * dist_value);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {1, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {2, 0.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(ex->eval(prov, {3, 0.0, 1.0}), 0.0);
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
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 2}, 2);
    tree.append(0, {0, 0, 10, 10}, {0, 0, 30, 5}, 3);

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
    tree.append(0, {0, 0, 10, 1}, {0, 0, 20, 2}, 2);
    tree.append(0, {0, 0, 10, 10}, {0, 0, 30, 5}, 3);

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

    // check operator
    e = thingify(2.0 + arb::iexpr::radius(3.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 + 3.0 * 10.0);

    // check unary operator
    e = thingify(+arb::iexpr::radius(3.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 3.0 * 10.0);
}

TEST(iexpr, sub) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::sub(arb::iexpr::scalar(2.0), arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 - 3.0 * 10.0);

    // check operator
    e = thingify(2.0 - arb::iexpr::radius(3.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 - 3.0 * 10.0);

    // check unary operator
    e = thingify(-arb::iexpr::radius(3.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), -3.0 * 10.0);
}

TEST(iexpr, mul) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::mul(arb::iexpr::scalar(2.0), arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 * 3.0 * 10.0);

    // check operator
    e = thingify(2.0 * arb::iexpr::radius(3.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 * 3.0 * 10.0);
}

TEST(iexpr, div) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::div(arb::iexpr::scalar(2.0), arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 / (3.0 * 10.0));

    // check operator
    e = thingify(2.0 / arb::iexpr::radius(3.0), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), 2.0 / (3.0 * 10.0));
}

TEST(iexpr, exp) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::exp(arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), std::exp(3.0 * 10.0));
}

TEST(iexpr, log) {
    segment_tree tree;
    tree.append(mnpos, {0, 0, 0, 10}, {0, 0, 10, 10}, 1);

    arb::mprovider prov(arb::morphology(std::move(tree)));

    auto e = thingify(arb::iexpr::log(arb::iexpr::radius(3.0)), prov);
    EXPECT_DOUBLE_EQ(e->eval(prov, {0, 0.0, 1.0}), std::log(3.0 * 10.0));
}
