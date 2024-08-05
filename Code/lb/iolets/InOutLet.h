// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLET_H
#define HEMELB_LB_IOLETS_INOUTLET_H

#include <memory>
#include "util/Vector3D.h"
#include "util/UnitConverter.h"
#include "lb/SimulationState.h"

namespace hemelb::redblood { class FlowExtension; }

namespace hemelb::lb
{
      //forward declare boundary comms class
      class BoundaryComms;
      class BoundaryCommunicator;

      /**
       * Base class for extra data needed by LB BC implementations.
       * Makes "Iolet coordinates" available.
       * These are coordinates in a frame aligned with the iolet plane.
       * iolet(0, 0, 0) corresponds to the iolet's position in lattice coordinates
       * x & y are arbitrary in plane components,
       * z is in the direction of the iolet's normal.
       */
      class InOutLet;
      class IoletExtraData
      {
        public:
          IoletExtraData(InOutLet const& iolet);
          virtual ~IoletExtraData();
          LatticePosition WorldToIolet(LatticeVector r);
          LatticePosition WorldToIolet(LatticePosition r);
        protected:
          using UnitVec = util::Vector3D<Dimensionless>;
          UnitVec e1;
          UnitVec e2;
          UnitVec n;
          LatticePosition centre;
        private:
          IoletExtraData();
      };

      /**
       * Base Iolet class
       * Contains information configured from the xml config file, and calculates a density near itself for use in LB calculation
       * Provides maximum and minimum range of densities/pressures for use by steering.
       */
      class InOutLet
      {
        public:
          InOutLet() :
              comms(nullptr), extraData(nullptr)
          {
          }

          virtual ~InOutLet() = default;

          /***
           * Copy the InOutLet.
           * @return Pointer to new IOLet.
           */
          [[nodiscard]] virtual InOutLet* clone() const = 0;

          /***
           * This is a castable? virtual method, which is perhaps an anti-pattern
           * We should potentially use dynamic cast checks instead.
           * @return true if any comms were done.
           */
          virtual bool IsCommsRequired() const
          {
            return false;
          }

          /***
           * This is a castable? virtual method, which is perhaps an anti-pattern
           * We should potentially use dynamic cast checks instead.
           * @return true if any comms were done.
           */
          virtual bool IsRegistrationRequired() const
          {
            return false;
          }
          void SetComms(BoundaryComms * boundaryComms)
          {
            comms = boundaryComms;
          }
          BoundaryComms * GetComms() const
          {
            return comms;
          }
          /***
           * Carry out communication necessary
           */
          virtual void DoComms(const BoundaryCommunicator& bcComms, const LatticeTimeStep timeStep);

          /***
           * Set up the Iolet.
           * @param units a UnitConverter instance.
           */
          virtual void Initialise(const util::UnitConverter* unitConverter)
          {
          }

          /***
           * Get the minimum density, in lattice units
           * @return minimum density, in lattice units
           */
          virtual LatticeDensity GetDensityMin() const = 0;

          /***
           * Get the maximum density, in lattice units
           * @return maximum density, in lattice units
           */
          virtual LatticeDensity GetDensityMax() const = 0;

          /***
           * Get the minimum pressure, in lattice units
           * @return
           */
          LatticePressure GetPressureMin() const
          {
            return GetDensityMin() * Cs2;
          }

          /***
           * Get the maximum pressure, in lattice units
           * @return
           */
          LatticePressure GetPressureMax() const
          {
            return GetDensityMax() * Cs2;
          }

          /// @todo: #632 This method must be moved to InOutletPressure
          virtual LatticeDensity GetDensity(LatticeTimeStep time_step) const = 0;

          /// @todo: #632 Is this method ever implemented not empty?
          virtual void Reset(SimulationState& state) = 0;

          const LatticePosition& GetPosition() const
          {
            return position;
          }

          void SetPosition(const LatticePosition& x)
          {
            position = x;
          }

          /**
           * Set the normal of the InOutlet
           * @param newNormal
           */
          void SetNormal(const util::Vector3D<Dimensionless>& newNormal)
          {
            normal = newNormal.GetNormalised();
          }

          const util::Vector3D<Dimensionless>& GetNormal() const
          {
            return normal;
          }

          /**
           * Set the minimum density throughout the simulation.
           * @param minSimDensity
           */
          void SetMinimumSimulationDensity(LatticeDensity minSimDensity)
          {
            minimumSimulationDensity = minSimDensity;
          }

          IoletExtraData* GetExtraData()
          {
            return extraData;
          }

          void SetExtraData(IoletExtraData* ed)
          {
            extraData = ed;
          }

          void SetFlowExtension(std::shared_ptr<redblood::FlowExtension> flowExt) {
            this->flowExtension = flowExt;
          }

          std::shared_ptr<redblood::FlowExtension> GetFlowExtension() const {
            return this->flowExtension;
          }

        protected:
          LatticeDensity minimumSimulationDensity;
          LatticePosition position;
          util::Vector3D<Dimensionless> normal;
          BoundaryComms* comms;
          IoletExtraData* extraData;
          std::shared_ptr<redblood::FlowExtension> flowExtension;
          friend class IoletExtraData;
      };

}

#endif /* HEMELB_LB_IOLETS_INOUTLET_H */
