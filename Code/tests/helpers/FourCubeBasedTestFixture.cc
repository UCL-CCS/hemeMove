// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/FourCubeBasedTestFixture.h"

#include "tests/helpers/OneInOneOutSimConfig.h"
#include "configuration/SimConfig.h"
#include "configuration/SimBuilder.h"

namespace hemelb::tests::helpers
{
    FourCubeBasedTestFixtureBase::FourCubeBasedTestFixtureBase(int cubesize)
            : initParams(), cubeSize(cubesize), cubeSizeWithHalo(cubesize + 2)
    {
        // +2 for the halo of empty valid locations around the cube
        latDat.reset(FourCubeLatticeData::Create(Comms(), cubesize + 2));
        dom = &latDat->GetDomain();
        OneInOneOutSimConfigReader reader;
        simConfig = reader.Read();
        simBuilder = std::make_unique<configuration::SimBuilder>(simConfig);

        simState = simBuilder->BuildSimulationState();
        lbmParams = lb::LbmParameters(simState->GetTimeStepLength(),
                                      simConfig.GetVoxelSize());
        unitConverter = simBuilder->GetUnitConverter();

        initParams.latDat = &latDat->GetDomain();
        initParams.lbmParams = &lbmParams;
        numSites = initParams.latDat->GetLocalFluidSiteCount();
    }

    FourCubeBasedTestFixtureBase::~FourCubeBasedTestFixtureBase() = default;

    lb::BoundaryValues FourCubeBasedTestFixtureBase::BuildIolets(
            geometry::SiteType tp,
            std::vector<configuration::IoletConfig> const& conf
    ) const {
        return {tp,
                *dom,
                simBuilder->BuildIolets(conf),
                simState.get(),
                Comms(),
                *unitConverter};
    }

    lb::BoundaryValues FourCubeBasedTestFixtureBase::BuildIolets(
            geometry::SiteType tp
    ) const {
        REQUIRE((tp == geometry::INLET_TYPE || tp == geometry::OUTLET_TYPE));
        auto conf = (tp == geometry::INLET_TYPE) ? simConfig.GetInlets() : simConfig.GetOutlets();
        return BuildIolets(tp, conf);
    }

}
