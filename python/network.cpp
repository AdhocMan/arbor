#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <arbor/version.hpp>
#include <arbor/arbexcept.hpp>
#include <arbor/common_types.hpp>
#include <arbor/network.hpp>

#include "strprintf.hpp"

namespace pyarb {
namespace py = pybind11;

void register_network(py::module& m) {
    using namespace py::literals;

    m.def("unique",
        &arb::unique,
        "s"_a,
        "Removes any duplicate global labels contained by merging overlapping intervals for "
        "matching local labels.");

    py::class_<arb::network_selection> network_selection(
        m, "network_selection", "Selects or rejects a connection when queried");

    network_selection
        .def_static("bernoulli_random",
            &arb::network_selection::bernoulli_random,
            "seed"_a,
            "p"_a,
            "Random selection using the bernoulli random distribution with probability \"p\" "
            "between 0.0 and 1.0")
        .def_static("custom",
            &arb::network_selection::custom,
            "func"_a,
            "Custom selection using the provided function \"func\". "
            "Repeated calls with the same arguments to \"func\" must yield the same result")
        .def_static("all", &arb::network_selection::all, "Select all")
        .def_static("none", &arb::network_selection::none, "Select none")
        .def_static("inter_cell",
            &arb::network_selection::all,
            "Only select connections between different cells")
        .def_static("not_equal",
            &arb::network_selection::not_equal,
            "Only select connections when the global labels are not equal. May select intra-cell "
            "connections, if the local label is not equal.")
        .def_static("invert", &arb::network_selection::invert, "s"_a, "Invert the selection")
        .def("__and__",
            [](const arb::network_selection& self, const arb::network_selection& other) {
                return self & other;
            })
        .def("__or__",
            [](const arb::network_selection& self, const arb::network_selection& other) {
                return self | other;
            })
        .def("__xor__",
            [](const arb::network_selection& self, const arb::network_selection& other) {
                return self ^ other;
            })
        .def(
            "__call__",
            [](const arb::network_selection& self,
                const arb::cell_global_label_type& src,
                const arb::cell_global_label_type& dest) { return self(src, dest); },
            "src"_a,
            "dest"_a);

    py::class_<arb::network_value> network_value(
        m, "network_value", "Provides a floating point value for a given connection");

    network_value.def(py::init([](double value) { return arb::network_value(value); }), "value"_a)
        .def_static("uniform_random",
            &arb::network_value::uniform_random,
            "seed"_a,
            "range"_a,
            "Uniform random value in (range[0], range[1]]. Always returns the same value for "
            "repeated  calls with the same arguments and calls are symmetric v(a, b) = v(b, a).")
        .def_static("custom",
            &arb::network_value::custom,
            "func"_a,
            "Custom value using the provided function \"func\". Repeated calls with the same "
            "arguments to \"func\" must yield the same result. For gap junction values, "
            "\"func\" must be symmetric (func(a,b) = func(b,a)).")
        .def_static("uniform",
            &arb::network_value::uniform,
            "value"_a,
            "Uniform value. Will always return the same value given at construction.")
        .def(
            "__call__",
            [](const arb::network_value& self,
                const arb::cell_global_label_type& src,
                const arb::cell_global_label_type& dest) { return self(src, dest); },
            "src"_a,
            "dest"_a);

    py::class_<arb::cell_connection_network> cell_connection_network(
        m, "cell_connection_network", "Generate reproducible list of cell connections.");

    cell_connection_network
        .def(py::init([](arb::network_value weight,
                          arb::network_value delay,
                          arb::network_selection selection,
                          arb::network_population src_pop,
                          arb::network_population dest_pop) {
            return arb::cell_connection_network(std::move(weight),
                std::move(delay),
                std::move(selection),
                std::move(src_pop),
                std::move(dest_pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "src_pop"_a,
            "dest_pop"_a)
        .def(py::init([](double weight,
                          arb::network_value delay,
                          arb::network_selection selection,
                          arb::network_population src_pop,
                          arb::network_population dest_pop) {
            return arb::cell_connection_network(weight,
                std::move(delay),
                std::move(selection),
                std::move(src_pop),
                std::move(dest_pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "src_pop"_a,
            "dest_pop"_a)
        .def(py::init([](arb::network_value weight,
                          double delay,
                          arb::network_selection selection,
                          arb::network_population src_pop,
                          arb::network_population dest_pop) {
            return arb::cell_connection_network(std::move(weight),
                delay,
                std::move(selection),
                std::move(src_pop),
                std::move(dest_pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "src_pop"_a,
            "dest_pop"_a)
        .def(py::init([](double weight,
                          double delay,
                          arb::network_selection selection,
                          arb::network_population src_pop,
                          arb::network_population dest_pop) {
            return arb::cell_connection_network(
                weight, delay, std::move(selection), std::move(src_pop), std::move(dest_pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "src_pop"_a,
            "dest_pop"_a)
        .def("generate", &arb::cell_connection_network::generate, "gid"_a)
        .def_property_readonly("weight",
            &arb::cell_connection_network::weight,
            "The weight used for generated cell connections.")
        .def_property_readonly("delay",
            &arb::cell_connection_network::delay,
            "The delay used for generated cell connections.")
        .def_property_readonly("selection",
            &arb::cell_connection_network::selection,
            "The network selection used for generated cell connections.")
        .def_property_readonly("source_population",
            &arb::cell_connection_network::source_population,
            "The source population used for generated cell connections.")
        .def_property_readonly("destination_population",
            &arb::cell_connection_network::destination_population,
            "The destination_population used for generated cell connections.");

    py::class_<arb::gap_junction_network> gap_junction_network(
        m, "gap_junction_network", "Generate reproducible list of gap junctions.");

    gap_junction_network
        .def(py::init([](arb::network_value weight,
                          arb::network_selection selection,
                          arb::network_population src_pop,
                          arb::network_population dest_pop) {
            return arb::gap_junction_network(
                std::move(weight), std::move(selection), std::move(src_pop), std::move(dest_pop));
        }),
            "weight"_a,
            "selection"_a,
            "src_pop"_a,
            "dest_pop"_a)
        .def(py::init([](double weight,
                          arb::network_selection selection,
                          arb::network_population src_pop,
                          arb::network_population dest_pop) {
            return arb::gap_junction_network(
                weight, std::move(selection), std::move(src_pop), std::move(dest_pop));
        }),
            "weight"_a,
            "selection"_a,
            "src_pop"_a,
            "dest_pop"_a)
        .def("generate", &arb::gap_junction_network::generate, "gid"_a)
        .def_property_readonly("weight",
            &arb::gap_junction_network::weight,
            "The weight used for generated gap junctions.")
        .def_property_readonly("selection",
            &arb::gap_junction_network::selection,
            "The network selection used for generated gap junctions.")
        .def_property_readonly("source_population",
            &arb::gap_junction_network::source_population,
            "The source population used for generated gap junctions.")
        .def_property_readonly("destination_population",
            &arb::gap_junction_network::destination_population,
            "The destination_population used for generated gap junctions.");
}

}  // namespace pyarb
