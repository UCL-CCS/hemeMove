// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_INITIALCONDITION_H
#define HEMELB_LB_INITIALCONDITION_H

#include <optional>
#include <variant>

#include "geometry/FieldData.h"
#include "configuration/SimConfig.h"

namespace hemelb::lb {

    class SimulationState;

    // This base class is a friend of SimulationState so it can set
    // the timestep
    struct InitialConditionBase {
      InitialConditionBase();
      InitialConditionBase(std::optional<LatticeTimeStep> t);

      void SetTime(SimulationState* sim) const;

    protected:
      // Allow access for derived classes (this is a friend of domain_type)
      inline distribn_t* GetFOld(geometry::FieldData* ld, site_t i) const {
	return ld->GetFOld(i);
      }
      inline distribn_t* GetFNew(geometry::FieldData* ld, site_t i) const {
	return ld->GetFNew(i);
      }

      mutable std::optional<LatticeTimeStep> initial_time;
    };
    
    struct EquilibriumInitialCondition : InitialConditionBase {
      
      EquilibriumInitialCondition();
      
      EquilibriumInitialCondition(std::optional<LatticeTimeStep> t0,
				  distribn_t rho,
				  distribn_t mx = 0.0, distribn_t my = 0.0, distribn_t mz = 0.0);
      
      template<class LatticeType>
      void SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const;
      
    private:
      distribn_t density;
      distribn_t mom_x;
      distribn_t mom_y;
      distribn_t mom_z;
    };
    
    struct CheckpointInitialCondition : InitialConditionBase {
      CheckpointInitialCondition(std::optional<LatticeTimeStep> t0, std::filesystem::path cp, std::optional<std::filesystem::path> maybeOff);
      
      template<class LatticeType>
      void SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const;

    private:
      std::filesystem::path cpFile;
      std::optional<std::filesystem::path> maybeOffFile;
    };

    class InitialCondition : std::variant<EquilibriumInitialCondition, CheckpointInitialCondition> {
      // Alias for private base
      using ICVar = std::variant<EquilibriumInitialCondition, CheckpointInitialCondition>;
    public:
      // Expose c'tors to allow creation
      using ICVar::ICVar;

      // Factory function for InitialCondition
      static InitialCondition FromConfig(const configuration::ICConfig&);
      
      void SetTime(SimulationState* sim) const;
      template<class LatticeType>
      void SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const;
    private:
      
      //ICVar ic;
    };

}
#endif // HEMELB_LB_INITIALCONDITION_H
