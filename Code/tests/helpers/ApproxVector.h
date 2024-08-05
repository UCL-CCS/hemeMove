// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_APPROXVECTOR_H
#define HEMELB_TESTS_HELPERS_APPROXVECTOR_H

#include <concepts>
#include <type_traits>

#include "util/Vector3D.h"

namespace hemelb::tests {
    // Similar to Catch2's Approx class, but for HemeLB Vector3D.
    // 
    // Don't use this class directly, use the template alias
    // `ApproxVector` below.
    //
    // Can construct with any valid arguments for the equivalent Vector3D constructors
    //
    // Compares the same if |this - other| < margin (default == 1e-6)
    template <std::floating_point T>
    struct ApproxVectorImpl {
        using vec = util::Vector3D<T>;
        using value_type = T;
    private:
        vec value;
        double margin = 1e-6;
    public:
        template <typename... ArgTs>
        explicit ApproxVectorImpl(ArgTs... args) : value(std::forward<ArgTs>(args)...) {
        }

        ApproxVectorImpl& Margin(double m) {
            margin = m;
            return *this;
        }

        // Factory for same margin, different value.
        template <typename... ArgTs>
        ApproxVectorImpl operator()(ArgTs... args) const {
            return ApproxVectorImpl{std::forward<ArgTs>(args)...}.Margin(margin);
        }

        friend bool operator==(const ApproxVectorImpl& lhs, const vec& rhs) {
            return lhs.impl(rhs);
        }
        friend bool operator==(const vec& lhs, const ApproxVectorImpl& rhs) {
            return rhs.impl(lhs);
        }
      
        friend bool operator!=(const ApproxVectorImpl& lhs, const vec& rhs) {
            return !(lhs == rhs);
        }
        friend bool operator!=(const vec& lhs, const ApproxVectorImpl& rhs) {
            return !(lhs == rhs);
        }

        friend std::ostream& operator<<(std::ostream& o, ApproxVectorImpl const& av)
        {
            return o << av.value << "+/-" << av.margin;
        }

    private:
        bool impl(const vec& other) const {
            return (value - other).GetMagnitudeSquared() < margin * margin;
        }
    };
    
    namespace detail {
        // build an overload set
        template <std::floating_point T>
        ApproxVectorImpl<T> help(T);
        template <std::floating_point T>
        ApproxVectorImpl<T> help(util::Vector3D<T>&&);
    }

    // This alias template allows you to pass either the element type
    // (e.g. double) or the Vector3D type
    // (e.g. hemelb::util::Vector3D<double>) and get the correct
    // ApproxVectorImpl out at the end.
    template <typename T>
    using ApproxVector = decltype(detail::help(std::declval<T>()));

    template <std::floating_point T>
    ApproxVectorImpl<T> ApproxV(util::Vector3D<T> v) {
        return ApproxVectorImpl<T>{std::move(v)};
    }

    namespace selftest {
        static_assert(std::is_same<ApproxVector<double>, ApproxVectorImpl<double>>::value, "ApproxVector helpers types all wrong");
        static_assert(std::is_same<ApproxVector<util::Vector3D<double>>, ApproxVectorImpl<double>>::value, "ApproxVector helpers types all wrong");
    }
}
#endif
