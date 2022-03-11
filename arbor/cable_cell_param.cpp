#include <cfloat>
#include <cmath>
#include <memory>
#include <numeric>
#include <vector>
#include <variant>
#include <tuple>

#include <arbor/cable_cell.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/s_expr.hpp>

#include "util/maputil.hpp"
#include "iexpr_impl.hpp"

namespace arb {

void check_global_properties(const cable_cell_global_properties& G) {
    auto& param = G.default_parameters;

    if (!param.init_membrane_potential) {
        throw cable_cell_error("missing global default parameter value: init_membrane_potential");
    }

    if (!param.temperature_K) {
        throw cable_cell_error("missing global default parameter value: temperature");
    }

    if (!param.axial_resistivity) {
        throw cable_cell_error("missing global default parameter value: axial_resistivity");
    }

    if (!param.membrane_capacitance) {
        throw cable_cell_error("missing global default parameter value: membrane_capacitance");
    }

    for (const auto& ion: util::keys(G.ion_species)) {
        if (!param.ion_data.count(ion)) {
            throw cable_cell_error("missing ion defaults for ion "+ion);
        }
    }

    for (const auto& kv: param.ion_data) {
        auto& ion = kv.first;
        const cable_cell_ion_data& data = kv.second;
        if (!data.init_int_concentration) {
            throw cable_cell_error("missing init_int_concentration for ion "+ion);
        }
        if (!data.init_ext_concentration) {
            throw cable_cell_error("missing init_ext_concentration for ion "+ion);
        }
        if (!data.init_reversal_potential && !param.reversal_potential_method.count(ion)) {
            throw cable_cell_error("missing init_reversal_potential or reversal_potential_method for ion "+ion);
        }
    }
}

cable_cell_parameter_set neuron_parameter_defaults = {
    // initial membrane potential [mV]
    -65.0,
    // temperatue [K]
    6.3 + 273.15,
    // axial resistivity [Ω·cm]
    35.4,
    // membrane capacitance [F/m²]
    0.01,
    // ion defaults:
    // internal concentration [mM], external concentration [mM], reversal potential [mV]
    {
        {"na", {10.0,  140.0,  115 - 65.}},
        {"k",  {54.4,    2.5,  -12 - 65.}},
        {"ca", {5e-5,    2.0,  12.5*std::log(2.0/5e-5)}}
    },
};


std::vector<defaultable> cable_cell_parameter_set::serialize() const {
    std::vector<defaultable> D;
    if (init_membrane_potential) {
        D.push_back(arb::init_membrane_potential{*this->init_membrane_potential});
    }
    if (temperature_K) {
        D.push_back(arb::temperature_K{*this->temperature_K});
    }
    if (axial_resistivity) {
        D.push_back(arb::axial_resistivity{*this->axial_resistivity});
    }
    if (membrane_capacitance) {
        D.push_back(arb::membrane_capacitance{*this->membrane_capacitance});
    }

    for (const auto& [name, data]: ion_data) {
        if (data.init_int_concentration) {
            D.push_back(init_int_concentration{name, *data.init_int_concentration});
        }
        if (data.init_ext_concentration) {
            D.push_back(init_ext_concentration{name, *data.init_ext_concentration});
        }
        if (data.init_reversal_potential) {
            D.push_back(init_reversal_potential{name, *data.init_reversal_potential});
        }
    }

    for (const auto& [name, mech]: reversal_potential_method) {
        D.push_back(ion_reversal_potential_method{name, mech});
    }

    if (discretization) {
        D.push_back(*discretization);
    }

    return D;
}

void decor::paint(region where, paintable what) {
    paintings_.push_back({std::move(where), std::move(what)});
}

void decor::place(locset where, placeable what, cell_tag_type label) {
    placements_.push_back({std::move(where), std::move(what), std::move(label)});
}

void decor::set_default(defaultable what) {
    std::visit(
            [this] (auto&& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<init_membrane_potential, T>) {
                    defaults_.init_membrane_potential = p.value;
                }
                else if constexpr (std::is_same_v<axial_resistivity, T>) {
                    defaults_.axial_resistivity = p.value;
                }
                else if constexpr (std::is_same_v<temperature_K, T>) {
                    defaults_.temperature_K = p.value;
                }
                else if constexpr (std::is_same_v<membrane_capacitance, T>) {
                    defaults_.membrane_capacitance = p.value;
                }
                else if constexpr (std::is_same_v<init_int_concentration, T>) {
                    defaults_.ion_data[p.ion].init_int_concentration = p.value;
                }
                else if constexpr (std::is_same_v<init_ext_concentration, T>) {
                    defaults_.ion_data[p.ion].init_ext_concentration = p.value;
                }
                else if constexpr (std::is_same_v<init_reversal_potential, T>) {
                    defaults_.ion_data[p.ion].init_reversal_potential = p.value;
                }
                else if constexpr (std::is_same_v<ion_reversal_potential_method, T>) {
                    defaults_.reversal_potential_method[p.ion] = p.method;
                }
                else if constexpr (std::is_same_v<cv_policy, T>) {
                    defaults_.discretization = std::forward<cv_policy>(p);
                }
            },
            what);
}

iexpr iexpr::scalar(double value) {
    return iexpr(iexpr_type::scalar, std::make_tuple(value));
}

iexpr iexpr::distance(double scale, mlocation loc) {
    return iexpr(iexpr_type::distance, std::make_tuple(scale, std::move(loc)));
}

iexpr iexpr::radius(double scale) {
    return iexpr(iexpr_type::radius, std::make_tuple(scale));
}

iexpr iexpr::diameter(double scale) {
    return iexpr(iexpr_type::diameter, std::make_tuple(scale));
}

iexpr iexpr::add(iexpr left, iexpr right) {
    return iexpr(iexpr_type::add, std::make_tuple(std::move(left), std::move(right)));
}

iexpr iexpr::mul(iexpr left, iexpr right) {
    return iexpr(iexpr_type::mul, std::make_tuple(std::move(left), std::move(right)));
}

std::unique_ptr<iexpr_interface> thingify(const iexpr& expr, const mprovider& m) {
    switch (expr.type_) {
        case iexpr_type::scalar:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::scalar(
                        std::get<0>(std::any_cast<std::tuple<double>>(expr))));
        case iexpr_type::distance:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::distance(
                        std::make_from_tuple<iexpr_impl::distance, std::tuple<double, mlocation>>(
                        std::any_cast<std::tuple<double, mlocation>>(expr))));
        case iexpr_type::radius:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::radius(
                        std::get<0>(std::any_cast<std::tuple<double>>(expr))));
        case iexpr_type::diameter:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::radius(
                        2 * std::get<0>(std::any_cast<std::tuple<double>>(expr))));
        case iexpr_type::add:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::add(
                        thingify(std::get<0>(std::any_cast<std::tuple<iexpr, iexpr>>(expr)), m),
                        thingify(std::get<1>(std::any_cast<std::tuple<iexpr, iexpr>>(expr)), m)));
        case iexpr_type::mul:
            return std::unique_ptr<iexpr_interface>(new iexpr_impl::add(
                        thingify(std::get<0>(std::any_cast<std::tuple<iexpr, iexpr>>(expr)), m),
                        thingify(std::get<1>(std::any_cast<std::tuple<iexpr, iexpr>>(expr)), m)));
    }
    return nullptr;
}

density::density(scaled_property<density> dens, const mprovider& provider) : density(std::move(dens.prop.mech)) {
    for(const auto& expr : dens.scale_expr) {
        scale_expr[expr.first] = thingify(expr.second, provider);
    }
}

} // namespace arb
