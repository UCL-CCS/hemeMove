// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REDBLOOD_PARALLEL_PARALLELFIXTURES_H
#define HEMELB_TESTS_REDBLOOD_PARALLEL_PARALLELFIXTURES_H

#include "redblood/Cell.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/FolderTestFixture.h"
#include "SimulationController.h"
#include "configuration/SimBuilder.h"

namespace hemelb::tests
{
    //! \brief gathers mid-domain and egde positions from all procs
    //! \details If there are insufficient number of edges, mid-domains are used instead.
    //! erase removes the components from the first process.
    std::vector<LatticeVector> GatherSpecialPositions(geometry::Domain const & domain,
						      size_t mid, size_t edges,
						      net::MpiCommunicator const &c);

    //! Make some functionality available
    class OpenedSimulationController : public SimulationController
    {
    public:
        //using SimulationController::SimulationController;
        using SimulationController::Finalise;
        OpenedSimulationController(SimulationController&& src) : SimulationController(src) {
        }

        void DoTimeStep()
        {
            SimulationController::DoTimeStep();
        }

        configuration::SimConfig const& GetSimConfig()
        {
            return this->simConfig;
        }
    };

    class OpenSimFixture : public helpers::FolderTestFixture {
    protected:
        template <typename STENCIL>
        using MyTraits = Traits<
                lb::DefaultLattice, lb::GuoForcingLBGK, lb::Normal,
                lb::DefaultStreamer, lb::DefaultWallStreamer, lb::DefaultInletStreamer, lb::DefaultOutletStreamer,
                STENCIL
        >;

        std::shared_ptr<configuration::CommandLine> options;
        template<class STENCIL>
        [[nodiscard]] auto CreateSim(net::IOCommunicator const &comm) const {
            using T = MyTraits<STENCIL>;
            auto tmp =  configuration::SimBuilder::CreateSim<T>(
                    *options, comm
            );
            return std::unique_ptr<OpenedSimulationController>(new OpenedSimulationController(std::move(*tmp)));
        }
    };

    class DummyCell : public NodeCell {
    public:
      LatticeForce force;

      DummyCell(LatticePosition const&position, LatticeForce f = 0e0,
		std::string const &templateName = "nope");
      DummyCell(std::vector<LatticePosition> const &positions, LatticeForce f = 0e0,
		std::string const &templateName = "nope");
      template<class ITER>
      DummyCell(ITER first, ITER last, LatticeForce f = 0e0, std::string const &templateName =
		"nope") :
	NodeCell(first, last, templateName), force(f)
      {
      }

      LatticeEnergy operator()(std::vector<LatticeForceVector> &f) const override;
      std::unique_ptr<CellBase> cloneImpl() const override;
    };

    //! Creates list of cells for each set of positions from each process
    std::vector<hemelb::redblood::CellContainer::value_type> CreateCellsFromSpecialPositions(
											     geometry::Domain const & domain, size_t mid, size_t edges,
											     net::MpiCommunicator const &c, size_t nCells);

    net::MpiCommunicator CreateDumbGraphComm(net::MpiCommunicator const &comm);
}

#endif  // ONCE
