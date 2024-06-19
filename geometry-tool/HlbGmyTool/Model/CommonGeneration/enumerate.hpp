// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_COMMON_ENUMERATE_HPP
#define HLBGMYTOOL_COMMON_ENUMERATE_HPP

namespace hemelb::gmytool {

template <typename Container>
class EnumIter {
 public:
  using iterator = decltype(std::begin(std::declval<Container>()));
  using iter_ref = typename std::iterator_traits<iterator>::reference;

 private:
  iterator iter_;
  int index_ = 0;

 public:
  explicit EnumIter(iterator begin) : iter_{begin}, index_{0} {}

  EnumIter& operator++() {
    ++iter_;
    ++index_;
    return *this;
  }

  friend bool operator==(const EnumIter& lhs, const EnumIter& rhs) {
    return lhs.iter_ == rhs.iter_;  // or self.index_ != rhs.index_;
  }
  friend bool operator!=(const EnumIter& lhs, const EnumIter& rhs) {
    return !(lhs == rhs);
  }

  std::pair<int, iter_ref> operator*() const { return {index_, *iter_}; }
};

template <typename Container>
class EnumerationAdaptor {
  using wrapper = EnumIter<Container>;

 public:
  EnumerationAdaptor(Container& container) : container_(container) {}
  wrapper begin() const { return wrapper{std::begin(container_)}; }
  wrapper end() const { return wrapper{std::end(container_)}; }

 private:
  Container& container_;
};

template <typename Container>
EnumerationAdaptor<Container> enumerate(Container& container) {
  return container;
}

template <typename Container>
EnumerationAdaptor<const Container> const_enumerate(
    const Container& container) {
  return container;
}

}  // namespace hemelb::gmytool
#endif
