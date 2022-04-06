#pragma once

#include <memory>
#include <optional>
#include <unordered_map>

#include <arbor/morph/locset.hpp>
#include <arbor/morph/region.hpp>
#include <arbor/iexpr.hpp>

namespace arb {

class label_dict {
    using ps_map = std::unordered_map<std::string, arb::locset>;
    using reg_map = std::unordered_map<std::string, arb::region>;
    using iexpr_map = std::unordered_map<std::string, arb::iexpr>;

    ps_map locsets_;
    reg_map regions_;
    iexpr_map iexpressions_;

public:
    void import(const label_dict& other, const std::string& prefix = "");

    void set(const std::string& name, locset ls);
    void set(const std::string& name, region reg);
    void set(const std::string& name, iexpr e);

    std::optional<arb::region> region(const std::string& name) const;
    std::optional<arb::locset> locset(const std::string& name) const;
    std::optional<arb::iexpr> iexpr(const std::string& name) const;

    const ps_map& locsets() const;
    const reg_map& regions() const;
    const iexpr_map& iexpressions() const;

    std::size_t size() const;
};

} //namespace arb
