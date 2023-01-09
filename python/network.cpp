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

    py::class_<arb::network_cell_group> network_cell_group(m,
        "network_cell_group",
        "For global identification of a range items.\n\n"
        "cell_global_range_label members:\n"
        "(1) the first global identifier of the range of cells [begin, end).\n"
        "(2) the past-the-end global identifier of the range of cells [begin, end).\n"
        "(3) a cell_local_label, referring to a labeled group of items on the cell and a policy "
        "for selecting a single item out of the group.\n");

    network_cell_group
        .def(py::init([](arb::cell_gid_type gid_begin,
                          arb::cell_gid_type gid_end,
                          std::vector<arb::cell_local_label_type> src_labels,
                          std::vector<arb::cell_local_label_type> dest_labels) {
            return arb::network_cell_group{
                gid_begin, gid_end, std::move(src_labels), std::move(dest_labels)};
        }),
            "gid_begin"_a,
            "gid_end"_a,
            "src_labels"_a,
            "dest_labels"_a,
            "Construct a network_cell_group identifier with [begin, end) global indices "
            "and a label argument identifying a item on a range of cells.\n"
            "The default round_robin policy is used for selecting one of possibly multiple items "
            "on the cell associated with the label.")
        .def(py::init([](py::tuple t) {
            if (py::len(t) != 4) throw std::runtime_error("tuple length != 4");
            return arb::network_cell_group{t[0].cast<arb::cell_gid_type>(),
                t[1].cast<arb::cell_gid_type>(),
                t[2].cast<std::vector<arb::cell_local_label_type>>(), t[3].cast<std::vector<arb::cell_local_label_type>>()};
        }),
            "Construct a cell_global_label identifier with tuple argument (begin, end, label):\n"
            "  begin:   The first global identifier of the range of cells [begin, end).\n"
            "  end:   The past-the-end global identifier of the range of cells [begin, end).\n"
            "  label: The cell_local_label representing the label and selection policy of an item "
            "on the cell.\n")
        .def_readwrite("gid_begin",
            &arb::network_cell_group::gid_begin,
            "The first global identifier of the range of cells [begin, end).")
        .def_readwrite("gid_end",
            &arb::network_cell_group::gid_end,
            "The past-the-end global identifier of the range of cells [begin, end).")
        .def_readwrite("src_labels",
            &arb::network_cell_group::src_labels,
            "The cell_local_label representing the label and selection policy of an item on the "
            "cell.")
        .def_readwrite("dest_labels",
            &arb::network_cell_group::dest_labels,
            "The cell_local_label representing the label and selection policy of an item on the "
            "cell.");
    // .def("__str__",
    //     [](arb::network_cell_group m) {
    //         return pprintf("<arbor.network_cell_group: begin {}, end {}, label ({}, {})>",
    //             m.begin,
    //             m.end,
    //             m.label.tag,
    //             m.label.policy);
    //     })
    // .def("__repr__", [](arb::network_cell_group m) {
    //     return pprintf("<arbor.network_cell_group: begin {}, end {}, label ({}, {})>",
    //         m.begin,
    //         m.end,
    //         m.label.tag,
    //         m.label.policy);
    // });

    py::class_<arb::spatial_network_cell_group> spatial_network_cell_group(m,
        "spatial_network_cell_group",
        "For global identification of a range items.\n\n"
        "cell_global_range_label members:\n"
        "(1) the first global identifier of the range of cells [begin, end).\n"
        "(2) the past-the-end global identifier of the range of cells [begin, end).\n"
        "(3) a cell_local_label, referring to a labeled group of items on the cell and a policy "
        "for selecting a single item out of the group.\n");

    spatial_network_cell_group
        .def(py::init([](arb::cell_gid_type gid_begin,
                          std::vector<arb::cell_local_label_type> src_labels,
                          std::vector<arb::cell_local_label_type> dest_labels,
                          std::vector<arb::network_location> locations) {
            return arb::spatial_network_cell_group{gid_begin,
                std::move(src_labels),
                std::move(dest_labels),
                std::move(locations)};
        }),
            "gid_begin"_a,
            "src_labels"_a,
            "dest_labels"_a,
            "locations"_a,
            "Construct a network_cell_group identifier with [begin, end) global indices "
            "and a label argument identifying a item on a range of cells.\n"
            "The default round_robin policy is used for selecting one of possibly multiple items "
            "on the cell associated with the label.")
        .def(py::init([](py::tuple t) {
            if (py::len(t) != 4) throw std::runtime_error("tuple length != 4");
            return arb::spatial_network_cell_group{t[0].cast<arb::cell_gid_type>(),
                t[1].cast<std::vector<arb::cell_local_label_type>>(),
                t[2].cast<std::vector<arb::cell_local_label_type>>(),
                t[3].cast<std::vector<arb::network_location>>()};
        }),
            "Construct a cell_global_label identifier with tuple argument (begin, end, label):\n"
            "  begin:   The first global identifier of the range of cells [begin, end).\n"
            "  end:   The past-the-end global identifier of the range of cells [begin, end).\n"
            "  label: The cell_local_label representing the label and selection policy of an item "
            "on the cell.\n");
        // .def_readwrite("gid_begin",
        //     &arb::network_cell_group::gid_begin,
        //     "The first global identifier of the range of cells [begin, end).")
        // .def_readwrite("gid_end",
        //     &arb::network_cell_group::gid_end,
        //     "The past-the-end global identifier of the range of cells [begin, end).")
        // .def_readwrite("src_labels",
        //     &arb::network_cell_group::src_labels,
        //     "The cell_local_label representing the label and selection policy of an item on the "
        //     "cell.")
        // .def_readwrite("dest_labels",
        //     &arb::network_cell_group::dest_labels,
        //     "The cell_local_label representing the label and selection policy of an item on the "
        //     "cell.");

    py::class_<arb::network_selection> network_selection(
        m, "network_selection", "Selects or rejects a connection when queried");

    py::class_<arb::spatial_network_selection> spatial_network_selection(
        m, "spatial_network_selection", "Selects or rejects a connection when queried");

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
            &arb::network_selection::inter_cell,
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
        .def("__and__",
            [](const arb::network_selection& self, const arb::spatial_network_selection& other) {
                return self & other;
            })
        .def("__or__",
            [](const arb::network_selection& self, const arb::network_selection& other) {
                return self | other;
            })
        .def("__or__",
            [](const arb::network_selection& self, const arb::spatial_network_selection& other) {
                return self | other;
            })
        .def("__xor__",
            [](const arb::network_selection& self, const arb::network_selection& other) {
                return self ^ other;
            })
        .def("__xor__",
            [](const arb::network_selection& self, const arb::spatial_network_selection& other) {
                return self ^ other;
            })
        .def(
            "__call__",
            [](const arb::network_selection& self,
                const arb::cell_global_label_type& src,
                const arb::cell_global_label_type& dest) { return self(src, dest); },
            "src"_a,
            "dest"_a);

    spatial_network_selection
        .def_static("custom",
            &arb::spatial_network_selection::custom,
            "func"_a,
            "Custom selection using the provided function \"func\". "
            "Repeated calls with the same arguments to \"func\" must yield the same result")
        .def_static("within_distance", &arb::spatial_network_selection::within_distance, "d"_a, "Select only within givin distance.")
        .def("__and__",
            [](const arb::spatial_network_selection& self,
                const arb::spatial_network_selection& other) { return self & other; })
        .def("__and__",
            [](const arb::spatial_network_selection& self, const arb::network_selection& other) {
                return self & other;
            })
        .def("__or__",
            [](const arb::spatial_network_selection& self,
                const arb::spatial_network_selection& other) { return self | other; })
        .def("__or__",
            [](const arb::spatial_network_selection& self, const arb::network_selection& other) {
                return self | other;
            })
        .def("__xor__",
            [](const arb::spatial_network_selection& self,
                const arb::spatial_network_selection& other) { return self ^ other; })
        .def("__xor__",
            [](const arb::spatial_network_selection& self, const arb::network_selection& other) {
                return self ^ other;
            })
        .def(
            "__call__",
            [](const arb::spatial_network_selection& self,
                const arb::cell_global_label_type& src,
                const arb::network_location& src_location,
                const arb::cell_global_label_type& dest,
                const arb::network_location& dest_location) {
                return self(src, src_location, dest, dest_location);
            },
            "src"_a,
            "src_location"_a,
            "dest"_a,
            "dest_location"_a);

    py::class_<arb::network_value> network_value(
        m, "network_value", "Provides a floating point value for a given connection");

    network_value.def(py::init([](double value) { return arb::network_value(value); }), "value"_a)
        .def_static("uniform_distribution",
            &arb::network_value::uniform_distribution,
            "seed"_a,
            "range"_a,
            "Uniform random value in (range[0], range[1]]. Always returns the same value for "
            "repeated  calls with the same arguments and calls are symmetric v(a, b) = v(b, a).")
        .def_static("normal_distribution",
            &arb::network_value::normal_distribution,
            "seed"_a,
            "mean"_a,
            "std_deviation"_a,
            "Radom value from a normal distribution with given mean and standard deviation. "
            "Always returns the same value for repeated calls with the same arguments and calls "
            "are symmetric v(a, b) = v(b, a).")
        .def_static("truncated_normal_distribution",
            &arb::network_value::truncated_normal_distribution,
            "seed"_a,
            "mean"_a,
            "std_deviation"_a,
            "range"_a,
            "Radom value from a truncated normal distribution with given mean and standard "
            "deviation (of a non-truncated normal distribution), where the value is always in "
            "(range[0], range[1]]. Always returns the same value for repeated calls with the same "
            "arguments and calls are symmetric v(a, b) = v(b, a). Note: Values are generated by "
            "reject-accept method from a normal distribution. Low acceptance rate can leed to poor "
            "performance, for example with very small ranges or a mean far outside the range.")
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
                          std::vector<arb::network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::spatial_network_value weight,
                          arb::spatial_network_value delay,
                          arb::spatial_network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::network_value weight,
                          arb::spatial_network_value delay,
                          arb::spatial_network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::network_value weight,
                          arb::network_value delay,
                          arb::spatial_network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::network_value weight,
                          arb::network_value delay,
                          arb::network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::network_value weight,
                          arb::spatial_network_value delay,
                          arb::network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::spatial_network_value weight,
                          arb::spatial_network_value delay,
                          arb::network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def(py::init([](arb::spatial_network_value weight,
                          arb::network_value delay,
                          arb::spatial_network_selection selection,
                          std::vector<arb::spatial_network_cell_group> pop) {
            return arb::cell_connection_network(
                std::move(weight), std::move(delay), std::move(selection), std::move(pop));
        }),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a)
        .def("generate", &arb::cell_connection_network::generate, "gid"_a);

    /*
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
    */
}

}  // namespace pyarb
