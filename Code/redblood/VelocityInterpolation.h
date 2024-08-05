// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_VELOCITYINTERPOLATION_H
#define HEMELB_REDBLOOD_VELOCITYINTERPOLATION_H

#include <iostream>

#include "units.h"
#include "Exception.h"
#include "geometry/FieldData.h"
#include "redblood/Interpolation.h"
#include "redblood/stencil.h"
#include "lb/kernels/GuoForcingLBGK.h"

#include <vector>

namespace hemelb::redblood
{
    //! \brief Returns the interpolated velocity at given point
    //! \param[in] center: off lattice position for which to interpolate.
    //! \param[in] latDat: lattice data
    //! \param[in] stencil: stencil to use for interpolation
    //! \tparam KERNEL: Needed to compute velocity
    template<class KERNEL>
    LatticeVelocity interpolateVelocity(geometry::FieldData const &latDat,
                                        LatticePosition const &center);

    namespace details
    {
      // Type traits to specialize for kernels with forces
      template<class KERNEL>
      struct HasForce : public std::false_type
      {
      };
      template<class LATTICE>
      struct HasForce<lb::GuoForcingLBGK<LATTICE> > : public std::true_type
      {
      };

      // Computes velocity for a given index on the lattice
      template<class KERNEL>
      struct VelocityFromLatticeData
      {
          VelocityFromLatticeData(geometry::FieldData const &latDat) :
              latticeData(latDat)
          {
          }
          // Computes velocity for given KERNEL and LatticeData
          LatticeVelocity operator()(LatticeVector const &indices) const
          {
            return operator()(latticeData.GetDomain().GetContiguousSiteId(indices));
          }
          LatticeVelocity operator()(size_t index) const
          {
            return computeVelocity(index, HasForce<KERNEL>());
          }

        protected:
          //! Implementation with forces correction
          LatticeVelocity computeVelocity(site_t index, std::true_type const &) const;
          //! Implementation without forces correction
          LatticeVelocity computeVelocity(site_t index, std::false_type const &) const;
          //! Lattice data that holds the grid of population and forces
          geometry::FieldData const &latticeData;
      };

      template<class KERNEL>
      LatticeVelocity VelocityFromLatticeData<KERNEL>::computeVelocity(
          site_t index, std::true_type const &) const
      {
        using LatticeType = typename KERNEL::LatticeType;
        LatticeVelocity result;
        distribn_t density;
        auto site = latticeData.GetSite(index);
        LatticeForceVector const &force(site.GetForce());
        auto const fDistribution = [&]() {
              if constexpr (build_info::USE_KRUEGER_ORDERING) {
                  // Use distribution functions at the beginning of the previous timestep (stored in
                  // FNew after the swap at the end of the timestep) in the IBM velocity interpolation.
                  // Follows approach in Timm's code
                  return latticeData.GetFNew<LatticeType>(index);
              } else {
                  return site.template GetFOld<LatticeType>();
              }
        }();
        LatticeType::CalculateDensityAndMomentum(fDistribution,
                                                 force,
                                                 density,
                                                 result);
        return result / density;
      }

      template<class KERNEL>
      LatticeVelocity VelocityFromLatticeData<KERNEL>::computeVelocity(
          site_t index, std::false_type const &) const
      {
        using LatticeType = typename KERNEL::LatticeType;
        LatticeMomentum mom;
        distribn_t density;
        auto const fDistribution = [&]() {
            if constexpr (build_info::USE_KRUEGER_ORDERING) {
                // Use distribution functions at the beginning of the previous timestep (stored in
                // FNew after the swap at the end of the timestep) in the IBM velocity interpolation.
                // Follows approach in Timm's code
                return latticeData.GetFNew<LatticeType>(index);
            } else {
                return latticeData.GetSite(index).template GetFOld<LatticeType>();
            }
        }();
        LatticeType::CalculateDensityAndMomentum(fDistribution,
                                                 density,
                                                 mom);
        return mom / density;
      }
    }

    template<class KERNEL, class STENCIL>
    LatticeVelocity interpolateVelocity(geometry::FieldData const &latDat,
                                        LatticePosition const &center)
    {
      auto iterator = interpolationIterator<STENCIL>(center);
      // Computes velocity for a given site index
      // Branches to one or another function depending on whether forces are available (since
      // velocity depends on forces)
      auto const gridfunc = details::VelocityFromLatticeData<KERNEL>(latDat);
      LatticeForceVector result(0, 0, 0);
      for (; iterator.IsValid(); ++iterator)
      {
        proc_t procid;
        site_t siteid;
        if (latDat.GetDomain().GetContiguousSiteId(*iterator, procid, siteid))
        {
          result += gridfunc(siteid) * iterator.weight();
        }
      }
      return result;
    }
} // hemelb::redblood
#endif
