#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <arbor/arbexcept.hpp>
#include <arbor/common_types.hpp>
#include <arbor/network.hpp>
#include <arbor/version.hpp>

#include "strprintf.hpp"

namespace pyarb {
namespace py = pybind11;

void register_network(py::module& m) {
    using namespace py::literals;

    py::class_<arb::network_cell_group> network_cell_group(m,
        "network_cell_group",
        "A group of cells and local labels identifying possible source and destination sites for "
        "cell connections.");

    network_cell_group
        .def(py::init([](arb::cell_gid_type gid_begin,
                          arb::cell_gid_type gid_end,
                          std::vector<arb::cell_local_label_type> src_labels,
                          std::vector<arb::cell_local_label_type> dest_labels,
                          std::vector<arb::cell_local_label_type> gj_labels) {
            return arb::network_cell_group{gid_begin,
                gid_end,
                std::move(src_labels),
                std::move(dest_labels),
                std::move(gj_labels)};
        }),
            "gid_begin"_a,
            "gid_end"_a,
            "src_labels"_a,
            "dest_labels"_a,
            "gj_labels"_a,
            "Construct a network_cell_group with [gid_begin, gid_end) global cell indices "
            "and a list of source and destination labels for cell connections, as well as a list "
            "of gap junction labels.\n")
        .def(py::init([](py::tuple t) {
            if (py::len(t) != 5) throw std::runtime_error("tuple length != 5");
            return arb::network_cell_group{t[0].cast<arb::cell_gid_type>(),
                t[1].cast<arb::cell_gid_type>(),
                t[2].cast<std::vector<arb::cell_local_label_type>>(),
                t[3].cast<std::vector<arb::cell_local_label_type>>(),
                t[4].cast<std::vector<arb::cell_local_label_type>>()};
        }),
            "Construct a network_cell_group with [gid_begin, gid_end) global cell indices "
            "and a list of source and destination labels, as well as a list of gap junction "
            "labels.\n")
        .def_readwrite("gid_begin",
            &arb::network_cell_group::gid_begin,
            "The first global identifier of the range of cells [begin, end).")
        .def_readwrite("gid_end",
            &arb::network_cell_group::gid_end,
            "The past-the-end global identifier of the range of cells [begin, end).")
        .def_readwrite("src_labels",
            &arb::network_cell_group::src_labels,
            "The source labels for cell connections.")
        .def_readwrite("dest_labels",
            &arb::network_cell_group::dest_labels,
            "The destination labels for cell connections.")
        .def_readwrite(
            "gj_labels", &arb::network_cell_group::dest_labels, "The gap junction labels.");

    py::class_<arb::spatial_network_cell_group> spatial_network_cell_group(m,
        "spatial_network_cell_group",
        "A group of cells and local labels identifying possible source and destination sites for "
        "cell connections, as well as a location for each cell.");

    spatial_network_cell_group
        .def(py::init([](arb::cell_gid_type gid_begin,
                          std::vector<arb::cell_local_label_type> src_labels,
                          std::vector<arb::cell_local_label_type> dest_labels,
                          std::vector<arb::cell_local_label_type> gj_labels,
                          std::vector<arb::network_location> locations) {
            return arb::spatial_network_cell_group{gid_begin,
                std::move(src_labels),
                std::move(dest_labels),
                std::move(gj_labels),
                std::move(locations)};
        }),
            "gid_begin"_a,
            "src_labels"_a,
            "dest_labels"_a,
            "gj_labels"_a,
            "locations"_a,
            "Construct a spatial_network_cell_group with global cell indices starting from "
            "gid_begin, a list of source / destination labels, a list of gap junction labels and a "
            "location for each cell.\n")
        .def(py::init([](py::tuple t) {
            if (py::len(t) != 5) throw std::runtime_error("tuple length != 5");
            return arb::spatial_network_cell_group{t[0].cast<arb::cell_gid_type>(),
                t[1].cast<std::vector<arb::cell_local_label_type>>(),
                t[2].cast<std::vector<arb::cell_local_label_type>>(),
                t[3].cast<std::vector<arb::cell_local_label_type>>(),
                t[4].cast<std::vector<arb::network_location>>()};
        }),
            "Construct a spatial_network_cell_group with global cell indices starting from "
            "gid_begin, a list of source / destination labels, a list of gap junction labels and a "
            "location for each cell.\n")
        .def_readonly("gid_begin",
            &arb::spatial_network_cell_group::gid_begin,
            "The first global identifier of the range of cells [begin, end).")
        .def_readonly("gid_end",
            &arb::spatial_network_cell_group::gid_end,
            "The past-the-end global identifier of the range of cells [begin, end).")
        .def_readwrite("src_labels",
            &arb::spatial_network_cell_group::src_labels,
            "The source labels for cell connections.")
        .def_readwrite("dest_labels",
            &arb::spatial_network_cell_group::dest_labels,
            "The destination labels for cell connections.")
        .def_readwrite(
            "gj_labels", &arb::spatial_network_cell_group::dest_labels, "The gap junction labels.")
        .def_readonly("locations",
            &arb::spatial_network_cell_group::locations,
            "The cell locations starting from gid_begin.");

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

    spatial_network_selection.def(py::init<arb::network_selection>())
        .def_static("custom",
            &arb::spatial_network_selection::custom,
            "func"_a,
            "Custom selection using the provided function \"func\". "
            "Repeated calls with the same arguments to \"func\" must yield the same result")
        .def_static("within_distance",
            &arb::spatial_network_selection::within_distance,
            "d"_a,
            "Select only within givin distance.")
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

    py::implicitly_convertible<arb::network_selection, arb::spatial_network_selection>();

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

    py::class_<arb::spatial_network_value> spatial_network_value(
        m, "spatial_network_value", "Provides a floating point value for a given connection");

    spatial_network_value.def(py::init<arb::network_value>());

    py::implicitly_convertible<arb::network_value, arb::spatial_network_value>();

    py::class_<arb::network_generator> network_generator(
        m, "network_generator", "Generate reproducible list of cell connections and gap junctions.");

    network_generator
        .def_static("cell_connections",
            static_cast<arb::network_generator (*)(arb::network_value,
                arb::network_value,
                arb::network_selection,
                std::vector<arb::network_cell_group>)>(arb::network_generator::cell_connections),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a,
            "Create a generator for cell connections.")
        .def_static("cell_connections",
            static_cast<arb::network_generator (*)(arb::spatial_network_value,
                arb::spatial_network_value,
                arb::spatial_network_selection,
                std::vector<arb::spatial_network_cell_group>)>(
                arb::network_generator::cell_connections),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "pop"_a,
            "Create a generator for spatial cell connections.")
        .def_static("gap_junctions",
            static_cast<arb::network_generator (*)(
                arb::network_value, arb::network_selection, std::vector<arb::network_cell_group>)>(
                arb::network_generator::gap_junctions),
            "gj_weight"_a,
            "gj_selection"_a,
            "pop"_a,
            "Create a generator for gap junctions.")
        .def_static("gap_junctions",
            static_cast<arb::network_generator (*)(arb::spatial_network_value,
                arb::spatial_network_selection,
                std::vector<arb::spatial_network_cell_group>)>(
                arb::network_generator::gap_junctions),
            "gj_weight"_a,
            "gj_selection"_a,
            "pop"_a,
            "Create a generator for gap junctions.")
        .def_static("combined",
            static_cast<arb::network_generator (*)(arb::network_value,
                arb::network_value,
                arb::network_selection,
                arb::network_value,
                arb::network_selection,
                std::vector<arb::network_cell_group>)>(arb::network_generator::combined),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "gj_weight"_a,
            "gj_selection"_a,
            "pop"_a,
            "Create a generator for cell connections and gap junctions.")
        .def_static("combined",
            static_cast<arb::network_generator (*)(arb::spatial_network_value,
                arb::spatial_network_value,
                arb::spatial_network_selection,
                arb::spatial_network_value,
                arb::spatial_network_selection,
                std::vector<arb::spatial_network_cell_group>)>(arb::network_generator::combined),
            "weight"_a,
            "delay"_a,
            "selection"_a,
            "gj_weight"_a,
            "gj_selection"_a,
            "pop"_a,
            "Create a generator for cell connections and gap junctions.");
}

}  // namespace pyarb
