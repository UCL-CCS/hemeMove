// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_QUANTITY_H
#define HEMELB_QUANTITY_H

#include <variant>

#include "Exception.h"
#include "util/ct_string.h"

namespace hemelb {

    // A simple wrapper around a representation type adding a units tag
    template<typename T, ct_string UNITS>
    struct quantity {
        T val;
    public:
        using type = T;
        static constexpr auto units = UNITS;

        quantity() = default;

        quantity(T x) : val(std::move(x)) {
        }

        T const &value() const {
            return val;
        }
    };

    // A variant of quantities having the same underlying type but different units.
    template<typename RepT, ct_string... possible_units>
    using quantity_union = std::variant<quantity<RepT, possible_units>...>;

    // Concept for the above
    template<typename T, typename VarT>
    concept IsVariantAlternative = requires(VarT v) {
        std::get<T>(v);
    };

    // Traits for the above.
    // Member `value` equivalent to being a quantity union.
    // representation_type is the common underlying value type.
    // N is the number of alternatives.
    // units is an array of the unit names.
    template<typename T>
    struct quantity_union_traits : std::false_type {
    };
    template<typename RepT, ct_string... possible_units>
    struct quantity_union_traits<std::variant<quantity<RepT, possible_units>...>> : std::true_type {
        using representation_type = RepT;
        static constexpr auto N = sizeof...(possible_units);
        static constexpr std::array<char const *, N> units = {possible_units.c_str()...};
    };

    template<typename T>
    concept QuantityUnion = quantity_union_traits<T>::value;

    // Helper to create quantity unions from a {value, units string} pair.
    // This is a bit tricky because of the interface between compile and run time code.

    // Primary template undefined
    template<typename T>
    struct quantity_union_factory;

    // Specialisation for actually being a quantity union
    template<typename RepT, ct_string... possible_units>
    struct quantity_union_factory<quantity_union<RepT, possible_units...>> {
        using type = std::variant<quantity<RepT, possible_units>...>;

        static constexpr auto N = sizeof...(possible_units);

        // Given the value and its units, construct a variant with
        // the appropriate active alternative or throw Exception.
        type operator()(RepT const &v, std::string_view units) {
            type ans;
            // usual trick for going from parameter pack to indices
            impl(std::make_index_sequence<N>{}, ans, v, units);
            return ans;
        }

        // Check every possible unit against the supplied units.
        // Standard trick of deducing indices with an index_sequence argument.
        // If none match then throw Exception.
        template<std::size_t... Is>
        void impl(std::index_sequence<Is...>, type &var, RepT const &v, std::string_view units) {
            // Do a fold over addition to get the sum of return values.
            // Each one returns N in the case of no match.
            auto sum_inds = (single<Is>(var, v, units) + ...);
            if (sum_inds == N * N) {
                throw Exception() << "No match for '" << units << "' in quantity union";
            }
        }

        // This checks if the given units match the units at index I
        // in the parameter pack. If yes, set the variant to be the
        // corresponding quantity, holding the given value.
        //
        // Return either I or N to distinguish the cases
        template<std::size_t I>
        auto single(type &var, RepT const &v, std::string_view units) {
            using this_quant = std::variant_alternative_t<I, type>;
            if (this_quant::units.c_str() == units) {
                var = this_quant{v};
                return I;
            } else {
                return N;
            }
        }

    };

    // Helper function for creating unions.
    template<typename RepT, ct_string... possible_units>
    auto make_quantity_union(RepT const &val, std::string_view units) {
        using U = quantity_union<RepT, possible_units...>;
        return quantity_union_factory<U>()(val, units);
    }
}
#endif