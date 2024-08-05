// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_BOUNDARYVALUES_H
#define HEMELB_LB_IOLETS_BOUNDARYVALUES_H

#include "net/IOCommunicator.h"
#include "net/IteratedAction.h"
#include "lb/iolets/InOutLet.h"
#include "lb/iolets/BoundaryCommunicator.h"
#include "geometry/SiteType.h"
#include "util/clone_ptr.h"

namespace hemelb::geometry { class Domain; }
namespace hemelb::lb
{

    class BoundaryValues : public net::IteratedAction
    {
        using IoletPtr = util::clone_ptr<InOutLet>;
    public:
        BoundaryValues(geometry::SiteType ioletType, geometry::Domain const& latticeData,
                       const std::vector<IoletPtr>& iolets,
                       SimulationState* simulationState, const net::MpiCommunicator& comms,
                       const util::UnitConverter& units);

        void RequestComms() override;
        void EndIteration() override;
        void Reset();

        void FinishReceive();

        LatticeDensity GetBoundaryDensity(const int index);

        LatticeDensity GetDensityMin(int boundaryId);
        LatticeDensity GetDensityMax(int boundaryId);

        static proc_t GetBCProcRank();

        // Borrow the pointer to an Iolet - this object still owns
        // the value.
        inline InOutLet* GetLocalIolet(unsigned int index)
        {
            return iolets[localIoletIDs[index]].get();
        }
        inline InOutLet const* GetLocalIolet(unsigned int index) const
        {
            return iolets[localIoletIDs[index]].get();
        }
        inline InOutLet const* GetGlobalIolet(unsigned int index) const
        {
            return iolets[index].get();
        }
        inline InOutLet* GetGlobalIolet(unsigned int index)
        {
            return iolets[index].get();
        }
        inline auto GetLocalIoletCount() const
        {
            return localIoletIDs.size();
        }
        inline auto GetGlobalIoletCount() const
        {
            return iolets.size();
        }
        inline unsigned int GetTimeStep() const
        {
            return state->GetTimeStep();
        }
        inline geometry::SiteType GetIoletType() const
        {
            return ioletType;
        }

    private:
        bool IsIoletOnThisProc(geometry::Domain const& latticeData, int boundaryId);
        std::vector<int> GatherProcList(bool hasBoundary);
        void HandleComms(InOutLet* iolet);
        geometry::SiteType ioletType;
        // All inlets/outlets in the simulation.
        // (Has to be a vector of pointers for InOutLet polymorphism)
        std::vector<IoletPtr> iolets;
        // The indices of the iolets that this process cares about
        std::vector<int> localIoletIDs;
        SimulationState* state;
        BoundaryCommunicator bcComms;
    };
}

#endif // HEMELB_LB_IOLETS_BOUNDARYVALUES_H
