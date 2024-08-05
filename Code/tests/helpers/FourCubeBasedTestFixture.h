// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H

#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

#include "configuration/SimConfig.h"
#include "lb/Collisions.h"
#include "lb/SimulationState.h"
#include "net/IOCommunicator.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/FolderTestFixture.h"


namespace hemelb
{
    namespace configuration {
        class SimBuilder;
    }
    namespace tests::helpers
    {

        class FourCubeBasedTestFixtureBase : public FolderTestFixture {

        public:
            FourCubeBasedTestFixtureBase(int cubesize);
            ~FourCubeBasedTestFixtureBase();

        protected:
            lb::BoundaryValues BuildIolets(geometry::SiteType tp) const;
            lb::BoundaryValues BuildIolets(geometry::SiteType tp, std::vector<configuration::IoletConfig> const& conf) const;
            std::unique_ptr<FourCubeLatticeData> latDat;
            // Non-owning ptr (owned by latDat)
            FourCubeDomain* dom;
            lb::InitParams initParams;
            site_t numSites;
            lb::LbmParameters lbmParams;
            configuration::SimConfig simConfig;
            std::unique_ptr<configuration::SimBuilder> simBuilder;
            std::shared_ptr<lb::SimulationState> simState;
            std::shared_ptr<const util::UnitConverter> unitConverter;
            int cubeSize;
            int cubeSizeWithHalo;
        };

        template <int CUBESIZE = 4>
        class FourCubeBasedTestFixture : public FourCubeBasedTestFixtureBase {
        public:
            FourCubeBasedTestFixture() : FourCubeBasedTestFixtureBase(CUBESIZE) {
            }
        };
    }
}
#endif
