// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H
#define HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H

#include "extraction/IterableDataSource.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/UnitConverter.h"

namespace hemelb::geometry { class FieldData; }
namespace hemelb::extraction
{
    class LbDataSourceIterator : public IterableDataSource
    {
    public:

        /**
         * Constructs an iterator using the LBM as a source.
         * @param propertyCache
         * @param data
         * @return
         */
        LbDataSourceIterator(const lb::MacroscopicPropertyCache& propertyCache,
                             const geometry::FieldData& data, int rank,
                             std::shared_ptr<util::UnitConverter> converter);

        /**
         * Reads the next fluid site from the data source. Returns true if values could
         * be obtained.
         *
         * @return
         */
        bool ReadNext() override;

        /**
         * Returns the coordinates of the site.
         * @return
         */
        [[nodiscard]] util::Vector3D<site_t> GetPosition() const override;

        /**
         * Returns the pressure at the site.
         * @return
         */
        [[nodiscard]] FloatingType GetPressure() const override;

        /**
         * Returns the velocity at the site.
         * @return
         */
        [[nodiscard]] util::Vector3D<FloatingType> GetVelocity() const override;

        /**
         * Returns the shear stress at the site.
         * @return
         */
        [[nodiscard]] FloatingType GetShearStress() const override;

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        [[nodiscard]] FloatingType GetVonMisesStress() const override;

        /**
         * Returns the shear rate at the site.
         * @return shear rate
         */
        [[nodiscard]] FloatingType GetShearRate() const override;

        /**
         * Returns the full stress tensor at the site.
         * @return stress tensor
         */
        [[nodiscard]] util::Matrix3D GetStressTensor() const override;

        /**
         * Returns the traction vector at a wall site (i.e. stress tensor times surface normal).
         * @return traction vector
         */
        [[nodiscard]] util::Vector3D<LatticeStress> GetTraction() const override;

        /**
         * Returns the projection of the traction vector on the tangential plane of a wall site.
         * @return projected traction vector
         */
        [[nodiscard]] util::Vector3D<LatticeStress> GetTangentialProjectionTraction() const override;

        /**
         * Returns a pointer to the velocity distribution of a site.
         * @return pointer to a velocity distribution
         */
        [[nodiscard]] const distribn_t* GetDistribution() const override;

        /**
         * Resets the iterator to the beginning again.
         */
        void Reset() override;

        /**
         * Returns true iff the passed location is within the lattice.
         *
         * @param
         * @return
         */
        [[nodiscard]] bool IsValidLatticeSite(const util::Vector3D<site_t>& location) const override;

        /**
         * Returns true iff the given location is available on this core (i.e. if the data
         * lives on this core).
         * @return
         */
        [[nodiscard]] bool IsAvailable(const util::Vector3D<site_t>& location) const override;

        /**
         * Returns the real-world size of a single lattice unit.
         * @return
         */
        [[nodiscard]] PhysicalDistance GetVoxelSize() const override;
        [[nodiscard]] PhysicalTime GetTimeStep() const override;
        [[nodiscard]] PhysicalMass GetMassScale() const override;
        [[nodiscard]] PhysicalPressure GetReferencePressure() const override;

        /**
         * Returns the origin of the geometry in real, spatial units.
         *
         * @return
         */
        [[nodiscard]] const PhysicalPosition& GetOrigin() const override;

        /**
         * Returns true if the site at the given location is marked as a wall site
         * (i.e. one of its links intersects a wall)
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        [[nodiscard]] bool IsWallSite(const util::Vector3D<site_t>& location) const override;

        /**
         * Returns the number of components in a velocity distribution.
         * @return
         */
        [[nodiscard]] unsigned GetNumVectors() const override;

      private:
        /**
         * The cache of properties for each site, which we iterate through.
         */
        const lb::MacroscopicPropertyCache& propertyCache;
        /**
         * The object containing information about the lattice.
         */
        const geometry::FieldData& data;
        /**
         * The rank of the current process in the LB communicator.
         */
        const int rank;
        /**
         * Object capable of converting from physical to lattice units and vice versa.
         */
        std::shared_ptr<util::UnitConverter> converter;
        /**
         * Iteration variable for tracking progress through all the local fluid sites.
         */
        site_t position;
    };
}

#endif /* HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H */
