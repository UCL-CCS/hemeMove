// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_VECTOR3D_H
#define HEMELB_UTIL_VECTOR3D_H

#include <array>
#include <cmath>
#include <iterator>
#include <ostream>
#include <algorithm>
#include <numeric>
#include <limits>
#include "util/numerical.h"

#ifdef HEMELB_CODE
#include "net/MpiDataType.h"
#endif

namespace hemelb::util
{
    /**
     * Base class for Vector3D that handles IndexError.
     *
     * Note that we DO NOT initialise the static member variable
     * HandlerFunction* handler in Vector3D.cc, as you might expect!
     * The reason is that we want to use different handlers within HemeLB
     * and the setup tool. For HemeLB proper, this initialisation is done in
     * Vector3DHemeLb.cc
     */
    class Vector3DBase
    {
      public:
        /**
         * Typedef for the error handler function.
         */
        using HandlerFunction = void (int);

        /**
         * Set the function used to handle errors, overriding the default.
         *
         * @param A pointer to the function to use to deal with an index
         * error. Must raise an exception or terminate the simulation or
         * similar.
         */
        static void SetIndexErrorHandler(HandlerFunction* func);
      protected:
        /**
         * Called by subclasses to handle an index error. Wraps the real
         * function used.
         *
         * @param The (out of range) index supplied
         */
        static HandlerFunction HandleIndexError;
      private:
        /**
         * The pointer to the real handler function to use.
         */
        static HandlerFunction* handler;
    };

    template <typename T>
    concept Arithmetic = std::integral<T> || std::floating_point<T>;

    // Helper traits
    // Note: standard defines additive (+, -) and multiplicative (*, /, %) categories of arithmetic operators
    //
    // The result type of doing T OPERATOR U
    template <Arithmetic T, Arithmetic U>
    using add_result_t = decltype(std::declval<std::decay_t<T>>() + std::declval<std::decay_t<U>>());
    template <Arithmetic T, Arithmetic U>
    using mul_result_t = decltype(std::declval<std::decay_t<T>>() * std::declval<std::decay_t<U>>());

    // Check for integer promotion shenanigans.
    //
    // This is used in the operators below to ensure you don't
    // accidentally have any surprises like how char + char -> int.
    // Fixes:
    // - explicitly cast one argument to the output type
    // - use a compound assignment (e.g +=)
    // - for Dot and Cross explicitly give the template args
    // - maybe remove the static assert and act like a builtin type?
    template <Arithmetic T, Arithmetic U, Arithmetic Res>
    constexpr bool integer_promotion_occurred = std::is_same_v<T, U> && !std::is_same_v<T, Res>;

    /**
     * Three dimensional vector class template.
     *
     * It exposes its members publicly for use and overloads arithmetic
     * operators as well as providing a number of useful methods. Other
     * methods are defined for convenience
     */
    template<Arithmetic T>
    class Vector3D : public Vector3DBase
    {
        // In debug mode, always do bounds checking.
        static constexpr bool debug =
#ifndef NDEBUG
                true
#else
                false
#endif
        ;

    public:
        // The type of the element
        using value_type = T;

        // The data.
        // Would like to make this private, but also want to allow
        // its use as a non-type template parameter, so this has to
        // be a structural type.
        std::array<value_type, 3> m_values;

        /**
         * Iterator pointing to first element of the Vector3D.
         * @return
         */
        auto begin()
        {
            return m_values.begin();
        }
        auto begin() const
        {
            return m_values.begin();
        }

        /**
         * Iterator pointing just past the last element of the Vector3D.
         * @return
         */
        auto end()
        {
            return m_values.end();
        }
        auto end() const
        {
            return m_values.end();
        }

        // Be friends with all other Vec3 instantiations.
        template<Arithmetic> friend class Vector3D;

        /**
         * x, y and z components.
         */
        constexpr T const& x() const {
            return m_values[0];
        }
        constexpr T & x() {
            return m_values[0];
        }
        constexpr T const& y() const {
            return m_values[1];
        }
        constexpr T & y() {
            return m_values[1];
        }
        constexpr T const& z() const {
            return m_values[2];
        }
        constexpr T & z() {
            return m_values[2];
        }

        // Default constructor. It's a value type so this has
        // undefined value when default constructed.
        constexpr Vector3D() = default;

        /**
         * Constructor accepting elements
         * @param x-component
         * @param y-component
         * @param z-component
         */
        constexpr Vector3D(const T iX, const T iY, const T iZ) :
            m_values({iX, iY, iZ})
        {
        }

        /**
         * Constructor filling in each component with the supplied argument.
         * @param used for all components
         */
        explicit constexpr Vector3D(const T iX) :
                m_values({iX, iX, iX})
        {
        }

        // Copy constructor from Vector3D<this type>
        constexpr Vector3D(const Vector3D&) = default;

        // Move constructor from Vector3D<this type>
        //
        // Note that this doesn't gain us anything, unless T is
        // usefully moveable. We also don't provide a converting move
        // c'tor as that will need a copy even if one of the types
        // involved is moveable.
        constexpr Vector3D(Vector3D&&) noexcept = default;

        // Constructor converting elements from another type
        template<class U>
        explicit constexpr Vector3D(const Vector3D<U>& i)
        {
            std::copy(i.m_values.begin(), i.m_values.end(), m_values.begin());
        }

        // Go the other way and create a copy as the supplied type
        template <Arithmetic U>
        constexpr Vector3D<U> as() const {
          return Vector3D<U>{*this};
        }

        // Copy and move assignment are both explicitly defaulted.
        // Deliberately do not want converting assign (use `as` instead).
        constexpr Vector3D& operator=(const Vector3D&) = default;
        constexpr Vector3D& operator=(Vector3D&&) noexcept = default;

        /**
         * Get a component by Direction with bounds check
         * @param lDirection
         * @return The component
         */
        constexpr T& at(std::size_t i)
        {
            try {
                return m_values.at(i);
            }
            catch (std::out_of_range& e) {
                // HandleIndexError will either throw or bring the simulation
                // down. We include this return statement to suppress warnings!
                HandleIndexError(i);
                return m_values[0];
            }
        }
        constexpr T const& at(std::size_t i) const
        {
            try {
                return m_values.at(i);
            }
            catch (std::out_of_range& e) {
                // HandleIndexError will either throw or bring the simulation
                // down. We include this return statement to suppress warnings!
                HandleIndexError(i);
                return m_values[0];
            }
        }
        /**
         * Get a component by index without bounds check in release builds
         * @param index
         * @return component
         */
        constexpr T& operator[](size_t index)
        {
            if constexpr (debug) {
                return at(index);
            } else {
                return m_values[index];
            }
        }

        // Getters through tuple interface to support structured binding.
        template <std::size_t I>
        constexpr auto get() {
            return at(I);
        }
        template <std::size_t I>
        constexpr auto get() const {
            return at(I);
        }

        /**
         * Get a component by index without bounds check in release builds
         * @param index
         * @return component
         */
        constexpr const T& operator[](size_t index) const
        {
            if constexpr (debug) {
                return at(index);
            } else {
                return m_values[index];
            }
        }

        /**
         * Equality
         */
        constexpr friend bool operator==(const Vector3D& left, const Vector3D& right)
        {
            return left.m_values == right.m_values;
        }
        constexpr friend bool operator!=(const Vector3D& left, const Vector3D& right)
        {
            return !(left == right);
        }

        /**
         * In-place normalisation (i.e. this will be a unit vector).
         * @return reference to this
         */
        Vector3D& Normalise()
        {
          (*this) /= GetMagnitude();
          return *this;
        }

        /**
         * Compute the unit vector that points in this direction.
         * @return the unit vector
         */
        constexpr Vector3D GetNormalised() const
        {
          Vector3D normed(*this);
          normed.Normalise();
          return normed;
        }

        /**
         * Compute the magnitude squared of the vector
         * @return magnitude**2
         */
        constexpr T GetMagnitudeSquared() const
        {
          return Dot(*this, *this);
        }

        /**
         * Compute the magnitude of the vector
         * @return magnitude
         */
        T GetMagnitude() const
        {
          static_assert(std::is_floating_point_v<T>,
                  "Vector3D::GetMagnitude only makes sense for floating types");
          return std::sqrt(GetMagnitudeSquared());
        }

        /**
         * In-place vector addition
         * @param right
         * @return the updated vector
         */
        template <Arithmetic U>
        constexpr Vector3D& operator+=(const Vector3D<U>& right)
        {
            for (std::size_t i = 0; i < m_values.size(); ++i)
                m_values[i] += right.m_values[i];
            return *this;
        }

        /**
         * Vector unary negation
         * @return -this
         */
        constexpr Vector3D operator-() const
        {
          return Vector3D(-x(), -y(), -z());
        }

        /**
         * In-place vector subtraction
         * @param right
         * @return this - right
         */
        template <Arithmetic U>
        constexpr Vector3D& operator-=(const Vector3D<U>& right)
        {
            for (std::size_t i = 0; i < m_values.size(); ++i)
                m_values[i] -= right.m_values[i];
            return *this;
        }

        /**
         * Multiplication by a scalar
         * @param multiplier
         * @return
         */
        template<Arithmetic U>
        friend constexpr auto operator*(const Vector3D& lhs, const U& rhs)
        {
            using RES = mul_result_t<T, U>;
            static_assert(!integer_promotion_occurred<T, U, RES>);
            return Vector3D<RES>{
                    lhs.x() * rhs,
                    lhs.y() * rhs,
                    lhs.z() * rhs
            };
        }

        // For scalar * vector just swap and call above
        template<Arithmetic U>
        friend constexpr auto operator*(const U& lhs, const Vector3D& rhs)
        {
            return rhs*lhs;
        }

        /**
         * In-place multiplication by a scalar
         * @param multiplier
         * @return
         */
        template<Arithmetic U>
        constexpr Vector3D& operator*=(const U& multiplier)
        {
            for (T& v: m_values)
                v *= multiplier;
            return *this;
        }

        // Division by a scalar
        // Return type is a new Vector3D of the type of (T / U)
        template<Arithmetic U>
        constexpr auto operator/(const U& divisor) const
        {
            using RES = mul_result_t<T, U>;
            static_assert(!integer_promotion_occurred<T, U, RES>);
            return Vector3D<RES>{x() / divisor,
                                 y() / divisor,
                                 z() / divisor};
        }

        // In-place division by a scalar
        // Returns the updated object
        template<Arithmetic U>
        constexpr Vector3D& operator/=(const U& divisor)
        {
            for (T& v: m_values)
                v /= divisor;
            return *this;
        }

        /**
         * Scalar modulus
         * @param divisor
         */
        template<Arithmetic U>
        constexpr auto operator%(const U& divisor) const
        {
            using RES = mul_result_t<T, U>;
            static_assert(!integer_promotion_occurred<T, U, RES>);
            return Vector3D<RES>{x() % divisor, y() % divisor, z() % divisor};
        }

        /**
         * In-place scalar modulus
         * @param divisor
         */
        template<Arithmetic U>
        constexpr Vector3D& operator%=(const U& divisor)
        {
            for (T& v: m_values)
                v %= divisor;
            return *this;
        }

        /**
         * Point-wise multiplication
         */
        constexpr Vector3D PointwiseMultiplication(const Vector3D& rightArgument) const
        {
          return {x() * rightArgument.x(), y() * rightArgument.y(), z() * rightArgument.z()};
        }

        /**
         * Point-wise division
         */
        constexpr Vector3D PointwiseDivision(const Vector3D& rightArgument) const
        {
          return {x() / rightArgument.x(), y() / rightArgument.y(), z() / rightArgument.z()};
        }

        /**
         * Updates the each component of this Vector3D with the smaller of the
         * corresponding component in this and the other Vector3D.
         *
         * @param vector to compare against
         */
        constexpr void UpdatePointwiseMin(const Vector3D& iCompareVector)
        {
            for (int i = 0; i < 3; ++i)
                m_values[i] = std::min(m_values[i], iCompareVector.m_values[i]);
        }

        /**
         * Updates the each component of this Vector3D with the larger of the
         * corresponding component in this and the other Vector3D.
         *
         * @param vector to compare against
         */
        constexpr void UpdatePointwiseMax(const Vector3D& iCompareVector)
        {
            for (int i = 0; i < 3; ++i)
                m_values[i] = std::max(m_values[i], iCompareVector.m_values[i]);
        }

        /**
         * Return whether the vector is within the bounds specified.
         * @param vMin
         * @param vMax
         * @return
         */
        constexpr bool IsInRange(const Vector3D& vMin, const Vector3D& vMax) const
        {
            bool ans = true;
            for (int i = 0; i < 3; ++i) {
                ans &= util::IsInRange(m_values[i], vMin.m_values[i], vMax.m_values[i]);
            }
            return ans;
        }

        /**
         * Vector filled with the maximum value for the element type.
         * @return
         */
        static constexpr Vector3D Largest()
        {
          return Vector3D(std::numeric_limits<T>::max());
        }

        /**
         * Vector filled with the minimum value for the element type.
         * @return
         */
        static constexpr Vector3D Lowest()
        {
          return Vector3D(std::numeric_limits<T>::lowest());
        }

        /**
         * Factory for Vector3Ds of ones.
         * @return
         */
        static constexpr Vector3D Ones()
        {
          return Vector3D(1);
        }

        /**
         * Factory for Vector3Ds of zeros.
         * @return
         */
        static constexpr Vector3D Zero()
        {
          return Vector3D(0);
        }
    };

    /**
      * Vector addition
      * @param right
      * @return this + right
     */
    template <Arithmetic T, Arithmetic U>
    constexpr auto operator+(const Vector3D<T>& left, const Vector3D<U>& right)
    {
        using RES = add_result_t<T, U>;
        static_assert(!integer_promotion_occurred<T, U, RES>);
        Vector3D<RES> ans;
        std::transform(left.begin(), left.end(), right.begin(), ans.begin(), std::plus<void>{});
        return ans;
    }

    /**
      * Vector subtraction
      * @param right
      * @return this - right
      */
    template <Arithmetic T, Arithmetic U>
    constexpr auto operator-(const Vector3D<T>& left, const Vector3D<U>& right)
    {
        using RES = add_result_t<T, U>;
        static_assert(!integer_promotion_occurred<T, U, RES>);
        Vector3D<RES> ans;
        std::transform(left.begin(), left.end(), right.begin(), ans.begin(), std::minus<void>{});
        return ans;
    }


    // Dot product
    template <Arithmetic T, Arithmetic U>
    constexpr auto Dot(Vector3D<T> const& l, Vector3D<U> const& r) {
        using RES = mul_result_t<T, U>;
        static_assert(!integer_promotion_occurred<T, U, RES>);
        return std::inner_product(l.m_values.begin(), l.m_values.end(), r.m_values.begin(), RES{0});
    }

    // Cross product
    template <Arithmetic T, Arithmetic U>
    constexpr auto Cross(Vector3D<T> const& l, Vector3D<U> const& r) {
        using RES = mul_result_t<T, U>;
        static_assert(!std::is_same_v<T, U> || std::is_same_v<T, RES>);
        return Vector3D<RES>{
                l.y() * r.z() - l.z() * r.y(),
                l.z() * r.x() - l.x() * r.z(),
                l.x() * r.y() - l.y() * r.x()
        };
    }

    template<Arithmetic T>
    std::ostream& operator<<(std::ostream& o, Vector3D<T> const& v3)
    {
      return o << "(" << v3.x() << "," << v3.y() << "," << v3.z() << ")";
    }

    namespace detail
    {
      bool CheckNextChar(std::istream& i, char c);
    }

    template<Arithmetic T>
    std::istream& operator>>(std::istream& i, Vector3D<T>& v3)
    {
      if (detail::CheckNextChar(i, '('))
      {
        i >> v3.x();
        if (detail::CheckNextChar(i, ','))
        {
          i >> v3.y();
          if (detail::CheckNextChar(i, ','))
          {
            i >> v3.z();
            detail::CheckNextChar(i, ')');
          }
        }
      }
      return i;
    }

#ifdef HEMELB_CODE
    // If we are in the main application, let MPI use the underlying array
    // for communications.
    template<typename T>
    MPI_Datatype MpiDataType(util::Vector3D<T> const& v) {
        return net::MpiDataType<std::array<T, 3>>();
    }
#endif
}

// Specialisations to support structured binding (officially sanctioned in the standard)
namespace std {
    template <class T>
    struct tuple_size< hemelb::util::Vector3D<T> >
      : std::integral_constant<std::size_t, 3>
    {};

    template <std::size_t I, class T>
    struct tuple_element<I, hemelb::util::Vector3D<T> >
    {
        using type = T;
    };
}

#endif // HEMELB_UTIL_VECTOR3D_H
