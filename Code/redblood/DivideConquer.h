//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_DIVIDE_AND_CONQUER_H
#define HEMELB_REDBLOOD_DIVIDE_AND_CONQUER_H

#include <vector>
#include <map>
#include <cmath>
#include "units.h"
#include "util/Vector3D.h"

namespace hemelb { namespace redblood {

namespace details { namespace {

  // Short-hand to get type of the base of Divide and Conquer class
  // Also aggregates key type and compare functor
  template<class T> struct DnCBase{
    typedef LatticeVector key_type;
    struct CompareKeys {
      bool operator()(key_type const &_a, key_type const &_b) const {
        if(_a.x > _b.x) return false;
        else if(_a.x < _b.x) return true;
        if(_a.y > _b.y) return false;
        else if(_a.y < _b.y) return true;
        return _a.z < _b.z;
      }
    };
    typedef std::multimap<key_type, T, CompareKeys> type;
  };

}}

//! \brief Multimap for divide and conquer algorithms
//! \details Items at a position x are mapped into boxes of a given size.
template<class T> class DivideConquer : public details::DnCBase<T>::type {
    typedef typename details::DnCBase<T>::type base_type;
  public:
    typedef typename base_type::key_type key_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef std::pair<iterator, iterator> range;
    typedef std::pair<const_iterator, const_iterator> const_range;

    //! Constructor sets size of cutoff
    DivideConquer(PhysicalDistance _boxsize) : base_type(), boxsize_(_boxsize) {}
    //! Insert into divide and conquer container
    iterator insert(LatticePosition const&_pos, T const &_value) {
      return base_type::insert(value_type(DowngradeKey(_pos), _value));
    }
    //! Insert into divide and conquer container
    iterator insert(key_type const&_pos, T const &_value) {
      return base_type::insert(value_type(_pos, _value));
    }
    //! All objects in a single divide and conquer box
    range equal_range(LatticePosition const &_pos) {
      return DivideConquer<T>::equal_range(DowngradeKey(_pos));
    }
    //! All objects in a single divide and conquer box
    range equal_range(key_type const &_pos) {
      return base_type::equal_range(_pos);
    }
    //! All objects in a single divide and conquer box
    const_range equal_range(LatticePosition const &_pos) const {
      return DivideConquer<T>::equal_range(DowngradeKey(_pos));
    }
    //! All objects in a single divide and conquer box
    const_range equal_range(key_type const &_pos) const {
      return base_type::equal_range(_pos);
    }

    //! Length of each box
    PhysicalDistance GetBoxSize() const { return boxsize_; }

    //! Converts from position to box index
    key_type DowngradeKey(LatticePosition const &_pos) const {
      return key_type(
        static_cast<LatticeCoordinate>(std::floor(_pos.x / boxsize_)),
        static_cast<LatticeCoordinate>(std::floor(_pos.y / boxsize_)),
        static_cast<LatticeCoordinate>(std::floor(_pos.z / boxsize_))
      );
    }
    //! No conversion since in box index type already
    key_type DowngradeKey(key_type const &_pos) const { return _pos; }
  protected:
    PhysicalDistance const boxsize_;
};

}} // hemelb::redblood

#endif
