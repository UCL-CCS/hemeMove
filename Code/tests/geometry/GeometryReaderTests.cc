// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>

#include "configuration/SimConfig.h"
#include "geometry/Domain.h"
#include "geometry/GeometryReader.h"
#include "lb/lattices/D3Q15.h"
#include "reporting/Timers.h"
#include "resources/Resource.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/LaddFail.h"
#include "tests/helpers/EqualitySiteData.h"

namespace hemelb::tests
{
    TEST_CASE_METHOD(helpers::FolderTestFixture, "GeometryReaderTests") {
      auto timings = std::make_unique<reporting::Timers>();
      
      auto reader = std::make_unique<geometry::GeometryReader>(lb::D3Q15::GetLatticeInfo(),
							       *timings,
							       Comms());
      CopyResourceToTempdir("four_cube.xml");
      CopyResourceToTempdir("four_cube.gmy");
      auto simConfig = configuration::SimConfig::New("four_cube.xml");

      SECTION("TestRead") {
	LADD_FAIL();
	reader->LoadAndDecompose(simConfig.GetDataFilePath());
      }

      SECTION("TestSameAsFourCube") {
	LADD_FAIL();
	auto fourCube = std::unique_ptr<FourCubeLatticeData>{FourCubeLatticeData::Create(Comms())};
	auto readResult = reader->LoadAndDecompose(simConfig.GetDataFilePath());
    auto&& dom = fourCube->GetDomain();

	for (site_t i = 1; i < 5; i++) {
	  for (site_t j = 1; j < 5; j++) {
	    bool isWallSite = (i == 1 || i == 4 || j == 1 || j == 4);

	    for (site_t k = 1; k < 5; k++) {
	      //std::cout << i << "," << j << "," << k << " > " << std::setbase(8) << fourCube->GetSiteData(i*16+j*4+k) << " : " << globalLattice->GetSiteData(i,j,k) << std::endl;
	      util::Vector3D<site_t> location{i, j, k};
	      site_t siteIndex = dom.GetGlobalNoncontiguousSiteIdFromGlobalCoords(location);

	      hemelb::geometry::SiteData siteData(readResult.Blocks[0].Sites[siteIndex]);
	      REQUIRE(fourCube->GetSite(dom.GetContiguousSiteId(location)).GetSiteData() == siteData);
	      //                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData().GetOtherRawData(),
	      //                                       siteData.GetOtherRawData());
	      //
	      //                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData().GetWallIntersectionData(),
	      //                                       siteData.GetWallIntersectionData());

	      REQUIRE(isWallSite == readResult.Blocks[0].Sites[siteIndex].wallNormalAvailable);

	      if (isWallSite) {
		/// @todo: #597 use CPPUNIT_ASSERT_EQUAL directly (having trouble with Vector3D templated over different types at the minute)
		/// CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetWallNormal(), readResult.Blocks[0].Sites[siteIndex].wallNormal);
		REQUIRE(fourCube->GetSite(dom.GetContiguousSiteId(location)).GetWallNormal()
			== readResult.Blocks[0].Sites[siteIndex].wallNormal.as<double>());
	      }
	    }
	  }
	}

      }

    }
}
