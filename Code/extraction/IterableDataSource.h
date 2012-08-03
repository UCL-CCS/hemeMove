// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_EXTRACTION_ITERABLEDATASOURCE_H
#define HEMELB_EXTRACTION_ITERABLEDATASOURCE_H

#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace extraction
  {
    class IterableDataSource
    {
      public:
        /**
         * Virtual destructor. We could declare this as pure virtual but that will cause a fail at
         * link-time (at least with GCC). GCC needs an object file to put the vtable in; by defining
         * this function in IterableDataSource.cc, we give it one.
         * @return
         */
        virtual ~IterableDataSource();

        /**
         * Reads the next fluid site from the data source. Its position,
         * pressure, velocity, etc are available from other functions.
         * Returns true if values could be obtained.
         *
         * @return
         */
        virtual bool ReadNext() = 0;

        /**
         * Returns the coordinates of the site.
         * @return
         */
        virtual util::Vector3D<site_t> GetPosition() const = 0;

        /**
         * Returns the pressure at the site.
         * @return
         */
        virtual float GetPressure() const = 0;

        /**
         * Returns the velocity at the site.
         * @return
         */
        virtual util::Vector3D<float> GetVelocity() const = 0;

        /**
         * Returns the shear stress at the site.
         * @return
         */
        virtual float GetShearStress() const = 0;

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        virtual float GetVonMisesStress() const = 0;

        /**
         * Returns the shear rate at the site.
         * @return
         */
        virtual float GetShearRate() const = 0;

        /**
         * Resets the iterator to the beginning again.
         */
        virtual void Reset() = 0;

        /**
         * Returns true iff the passed location is within the lattice.
         *
         * @param
         * @return
         */
        virtual bool IsValidLatticeSite(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns true iff the given location is available on this core (i.e. if the data
         * lives on this core).
         * @return
         */
        virtual bool IsAvailable(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns the real-world size of a single lattice unit.
         * @return
         */
        virtual distribn_t GetVoxelSize() const = 0;

        /**
         * Returns the origin of the geometry in real units.
         * @return
         */
        virtual const util::Vector3D<distribn_t>& GetOrigin() const = 0;

        /**
         * Returns true if the site at the given location is marked as an edge site
         * (i.e. one of its links intersects a wall)
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        virtual bool IsEdgeSite(const util::Vector3D<site_t>& location) const = 0;
    };
  }
}

#endif /* HEMELB_EXTRACTION_ITERABLEDATASOURCE_H */
