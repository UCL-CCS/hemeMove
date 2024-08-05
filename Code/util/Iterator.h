// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_ITERATOR_H
#define HEMELB_UTIL_ITERATOR_H

#include <iostream>
#include <tuple>
#include <compare>

namespace hemelb::util
{
    //! For-ranged loop plus enumeration
    template<std::forward_iterator WrappedIter, std::sentinel_for<WrappedIter> WrappedEnd, class CounterT = std::iter_difference_t<WrappedIter>>
    class Enumerate
    {
    protected:
        WrappedIter first;
        WrappedEnd last;
        using WrappedRefT = std::iter_reference_t<WrappedIter>;
        //! Holds info about enumeration
        struct EnumerateItem
        {
            CounterT index;
            WrappedRefT value;
        };

    public:
        struct sentinel {
            WrappedEnd sent;
        };

        class iterator
        {
        protected:
            CounterT index;
            WrappedIter iter;
        public:
            explicit iterator(WrappedIter const &iter) : index(0), iter(iter)
            {
            }
            iterator(WrappedIter const &iter, size_t index) : index(index), iter(iter)
            {
            }
            iterator& operator++()
            {
              ++index;
              ++iter;
              return *this;
            }
            iterator operator++(int) &
            {
              iterator result(iter.iter);
              iter++;
              return result;
            }
            EnumerateItem operator*() const
            {
              return {index, *iter};
            }

            friend bool operator==(const iterator& l, const iterator& r) {
                return l.iter == r.iter;
            }

            friend bool operator==(iterator const& it, sentinel const& s)
            requires (!std::same_as<WrappedIter, WrappedEnd>)
            {
                return it.iter == s.sent;
            }
        };

        Enumerate(WrappedIter const &first, WrappedEnd const &last) : first(first), last(last)
        {
        }
        iterator begin() const
        {
          return {first, 0};
        }
        auto end() const
        {
            if constexpr (std::same_as<WrappedIter, WrappedEnd>) {
                return iterator(last, std::distance(first, last));
            } else {
                return sentinel{last};
            }
        }
    };

    //! Ranged-for loops with enumeration
    template<class CONTAINER>
    auto enumerate(CONTAINER&& c)
    {
        using namespace std;
        auto b = begin(std::forward<CONTAINER>(c));
        auto e = end(std::forward<CONTAINER>(c));
        using IterT = decltype(b);
        using IterOrSentinelT = decltype(e);
        return Enumerate<IterT, IterOrSentinelT>(b, e);
    }

    //! Ranged-for loops with enumeration
    template<class CounterT, class CONTAINER>
    auto enumerate_with(CONTAINER&& c)
    {
        using namespace std;
        auto b = begin(std::forward<CONTAINER>(c));
        auto e = end(std::forward<CONTAINER>(c));
        using IterT = decltype(b);
        using IterOrSentinelT = decltype(e);
        return Enumerate<IterT, IterOrSentinelT, CounterT>(b, e);
    }

    //! Ranged-for loops with enumeration
    template<class CONTAINER>
    auto cenumerate(CONTAINER const& c)
    {
        using namespace std;
        auto b = begin(c);
        auto e = end(c);
        using IterT = decltype(b);
        using IterOrSentinelT = decltype(e);
        return Enumerate<IterT, IterOrSentinelT>(b, e);
    }

    //! For-ranged loop over several containers in parallel
    template<class... ContainerTs>
    class Zip
    {
    public:
        //! Iterates over all objects that are being zipped
        template<class ... ITERATOR>
        class zip_iterator
        {
            static constexpr auto N = sizeof...(ITERATOR);
        public:
            using reference = std::tuple<typename std::iterator_traits<ITERATOR>::reference...>;
            //! Inherits from all iterator types
            zip_iterator(ITERATOR const & ... args) : iter(args...)
            {
            }
            //! Inherits from all iterator types
            zip_iterator(ITERATOR && ... args) : iter(std::move(args)...)
            {
            }
            zip_iterator& operator++()
            {
                increment(std::make_index_sequence<N>());
                return *this;
            }
            reference operator*() const
            {
                return deref(std::make_index_sequence<N>());
            }

            //! True if any of the iterator is equal to the corresponding on in b
            friend bool operator==(zip_iterator const &a, zip_iterator const &b)
            {
                return IsEqual<0>(a, b);
            }
            friend bool operator!=(zip_iterator const& a, zip_iterator const &b)
            {
                return !(a == b);
            }

        protected:
            // For each of these helpers, we use an index sequence
            // as the argument to deduce the template arguments `Is`.
            // Should be called with `std::make_index_sequence<N>()`.

            //! Increment
            // Note folding over the comma operator
            template <std::size_t... Is>
            void increment(std::index_sequence<Is...>) {
                (++std::get<Is>(iter), ...);
            }
            //! Comparison
            // Do a fold over logical and
            template<std::size_t... Is>
            friend bool IsEqual(zip_iterator const& a, zip_iterator const &b)
            {
                return ((std::get<Is>(a.iter) == std::get<Is>(b.iter)) && ...);
            }
            //! Derefence
            template <std::size_t... Is>
            reference deref(std::index_sequence<Is...>) const {
                return {*std::get<Is>(iter)...};
            }

            //! Holds each iterator
            std::tuple<ITERATOR...> iter;
        };

        using iterator = zip_iterator<decltype(std::declval<ContainerTs>().begin())...>;

        Zip(ContainerTs && ... cont) : first(cont.begin()...), last(cont.end()...)
        {
        }
        iterator begin() const
        {
          return first;
        }
        iterator end() const
        {
          return last;
        }
      protected:
        iterator first;
        iterator last;
    };

    //! \brief Ranged-for loops over sets of containers in parallel
    //! \code
    //! for(auto const && item: zip(container0, container2, containerN))
    //! {
    //!   std::get<0>(item);
    //!   std::get<1>(item);
    //! }
    //! \endcode
    template<class ... CONTAINER>
    auto zip(CONTAINER && ... x) -> Zip<CONTAINER...>
    {
      return {x...};
    }
    //! \brief Ranged-for loops over sets of (const) containers in parallel
    //! \code
    //! for(auto const && item: zip(container0, container2, containerN))
    //! {
    //!   std::get<0>(item);
    //!   std::get<1>(item);
    //! }
    //! \endcode
    template<class ... CONTAINER>
    auto czip(CONTAINER const & ... x) -> Zip<CONTAINER const &...>
    {
      return {x...};
    }

}

#endif
