// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_ITERABLEDATASOURCE_H
#define HEMELB_EXTRACTION_ITERABLEDATASOURCE_H

#include "util/Vector3D.h"
#include "units.h"
#include "util/Matrix3D.h"

namespace hemelb::extraction
{

    using FloatingType = float;

    class IterableDataSource
    {
    public:
        /**
         * Virtual destructor. We could declare this as pure virtual but that will cause a fail at
         * link-time (at least with GCC). GCC needs an object file to put the vtable in; by defining
         * this function in IterableDataSource.cc, we give it one.
         * @return
         */
        virtual ~IterableDataSource() = default;

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
        [[nodiscard]] virtual util::Vector3D<site_t> GetPosition() const = 0;

        /**
         * Returns the pressure at the site.
         * @return
         */
        [[nodiscard]] virtual FloatingType GetPressure() const = 0;

        /**
         * Returns the velocity at the site.
         * @return
         */
        [[nodiscard]] virtual util::Vector3D<FloatingType> GetVelocity() const = 0;

        /**
         * Returns the shear stress at the site.
         * @return
         */
        [[nodiscard]] virtual FloatingType GetShearStress() const = 0;

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        [[nodiscard]] virtual FloatingType GetVonMisesStress() const = 0;

        /**
         * Returns the shear rate at the site.
         * @return
         */
        [[nodiscard]] virtual FloatingType GetShearRate() const = 0;

        /**
         * Returns the full stress tensor at the site.
         * @return stress tensor
         */
        [[nodiscard]] virtual util::Matrix3D GetStressTensor() const = 0;

        /**
         * Returns the traction vector at a wall site (i.e. stress tensor times surface normal).
         * @return traction vector
         */
        [[nodiscard]] virtual util::Vector3D<LatticeStress> GetTraction() const = 0;

        /**
         * Returns the projection of the traction vector on the tangential plane of a wall site.
         * @return projected traction vector
         */
        [[nodiscard]] virtual util::Vector3D<LatticeStress> GetTangentialProjectionTraction() const = 0;

        /**
         * Returns a pointer to the velocity distribution at a site.
         * @return pointer to velocity distribution
         */
        [[nodiscard]] virtual const distribn_t* GetDistribution() const = 0;

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
        [[nodiscard]] virtual bool IsValidLatticeSite(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns true iff the given location is available on this core (i.e. if the data
         * lives on this core).
         * @return
         */
        [[nodiscard]] virtual bool IsAvailable(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns the real-world size of a single lattice unit.
         * @return
         */
        [[nodiscard]] virtual PhysicalDistance GetVoxelSize() const = 0;
        [[nodiscard]] virtual PhysicalTime GetTimeStep() const = 0;
        [[nodiscard]] virtual PhysicalMass GetMassScale() const = 0;
        [[nodiscard]] virtual PhysicalPressure GetReferencePressure() const = 0;

        /**
         * Returns the origin of the geometry in real units.
         * @return
         */
        [[nodiscard]] virtual const PhysicalPosition& GetOrigin() const = 0;

        /**
         * Returns true if the site at the given location is marked as a wall site
         * (i.e. one of its links intersects a wall)
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        [[nodiscard]] virtual bool IsWallSite(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns the number of components in a velocity distribution.
         * @return
         */
        [[nodiscard]] virtual unsigned GetNumVectors() const = 0;
    };
}
#endif /* HEMELB_EXTRACTION_ITERABLEDATASOURCE_H */
