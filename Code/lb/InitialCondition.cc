// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/InitialCondition.h"

#include "lb/SimulationState.h"

namespace hemelb::lb {

    // InitialConditionBase
    
    InitialConditionBase::InitialConditionBase() {
    }
    InitialConditionBase::InitialConditionBase(std::optional<LatticeTimeStep> t) : initial_time(t) {
    }

    void InitialConditionBase::SetTime(SimulationState* sim) const {
        if (initial_time)
            sim->timeStep = *initial_time;
    }


    // Equilibrium
    
    EquilibriumInitialCondition::EquilibriumInitialCondition() :
      InitialConditionBase(),
      density(1.0), mom_x(0.0), mom_y(0.0), mom_z(0.0) {
    }
    
    EquilibriumInitialCondition::EquilibriumInitialCondition(
      std::optional<LatticeTimeStep> t0,
      distribn_t rho,
      distribn_t mx, distribn_t my, distribn_t mz) :
      InitialConditionBase(t0),
      density(rho),
      mom_x(mx), mom_y(my), mom_z(mz) {
    }
    
    CheckpointInitialCondition::CheckpointInitialCondition(std::optional<LatticeTimeStep> t0, std::filesystem::path cp, std::optional<std::filesystem::path> maybeOff)
      : InitialConditionBase(t0), cpFile(std::move(cp)), maybeOffFile(std::move(maybeOff)) {
    }

    // InitialCondition - sum type container

    // Visitor for setting time
    struct TSetter {
        using result_type = void;
        template <typename T>
        void operator()(T t) const {
            t.SetTime(ss);
        }
        SimulationState* ss;
    };
    void InitialCondition::SetTime(SimulationState* sim) const {
        const ICVar* self = this;
        std::visit(TSetter{sim}, *self);
    }

    // See InitialCondtions.hpp for setting Fs (distributions)
}
