// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <tinyxml2.h>

#include "configuration/SimConfigReader.h"
#include "redblood/FlowExtension.h"
#include "util/Vector3D.h"
#include "units.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/SimConfBuildHelp.h"

namespace hemelb::tests
{
    using namespace redblood;
    
    TEST_CASE_METHOD(FlowExtensionFixture, "FlowExtensionTests", "[redblood]") {
      using Point = util::Vector3D<LatticeDistance>;

      SECTION("testAxis") {
	// Check points along the axis from just outside the start of the cylinder (x = 0) to just inside
	INFO("Point (-0.001, 0, 0) is outside the cylinder");
	REQUIRE(!contains(flowExt, Point(-0.001, 0, 0)));
	INFO("Point ( 0.000, 0, 0) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.000, 0, 0)));
	INFO("Point ( 0.001, 0, 0) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.001, 0, 0)));

	// Check points along the axis from just inside the end of the cylinder (x = 10) to just outside
	INFO("Point ( 9.999, 0, 0) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(9.999, 0, 0)));
	INFO("Point (10.000, 0, 0) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(10.000, 0, 0)));
	INFO("Point (10.001, 0, 0) is outside the cylinder");
	REQUIRE(!contains(flowExt, Point(10.001, 0, 0)));
      }

      SECTION("testCircumference") {
	// Points around the circumference should be inside the cylinder
	INFO("Point (0.0,    1.000,   0.000 ) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, 1.0, 0.0)));
	INFO("Point (0.0,  cos(45),  sin(45)) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, cos(45), sin(45))));
	INFO("Point (0.0,    0.000,   1.000 ) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, 0.0, 1.0)));
	INFO("Point (0.0,  cos(45), -sin(45)) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, cos(45), -sin(45))));
	INFO("Point (0.0,   -1.000,   0.000 ) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, -1.0, 0.0)));
	INFO("Point (0.0, -cos(45), -sin(45)) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, -cos(45), -sin(45))));
	INFO("Point (0.0,    0.000,  -1.000 ) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, 0.0, -1.0)));
	INFO("Point (0.0, -cos(45),  sin(45)) is inside the cylinder");
	REQUIRE(contains(flowExt, Point(0.0, -cos(45), sin(45))));

	// Points forming a square around the axis should be outside the cylinder
	INFO("Point(0.0,  1.0,  1.0) is outside the cylinder");
	REQUIRE(!contains(flowExt, Point(0.0, 1.0, 1.0)));
	INFO("Point(0.0,  1.0, -1.0) is outside the cylinder");
	REQUIRE(!contains(flowExt, Point(0.0, 1.0, -1.0)));
	INFO("Point(0.0, -1.0, -1.0) is outside the cylinder");
	REQUIRE(!contains(flowExt, Point(0.0, -1.0, -1.0)));
	INFO("Point(0.0, -1.0,  1.0) is outside the cylinder");
	REQUIRE(!contains(flowExt, Point(0.0, -1.0, 1.0)));
      }

      auto approx = Approx(0.0).margin(1e-8);

      SECTION("testLinearWeight") {
	FlowExtension const flow(util::Vector3D<Dimensionless>(-1.0, 0, 0),
				 LatticePosition(0.5, 0.5, 0.5),
				 2.0,
				 0.5,
				 1.5);
	std::vector<Dimensionless> const ys{0.5, 0.7};
	std::vector<Dimensionless> const epsilons{0.3, 0.5, 0.7, 1.};
	for (auto& y : ys) {
	  REQUIRE(approx(0) == linearWeight(flow, LatticePosition(2.4, y, 0.5)));
	  REQUIRE(approx(1) == linearWeight(flow, LatticePosition(0.5, y, 0.5)));
	  for (auto& epsilon : epsilons) {
	    auto const pos = flow.origin + flow.normal * flow.fadeLength * epsilon
	      + LatticePosition(0, y - 0.5, 0);
	    REQUIRE(approx(1e0 - epsilon) == linearWeight(flow, pos));
	  }
	}
	REQUIRE(approx(0.0) == linearWeight(flow, LatticePosition(0.5, 0.5, 1.0 + 1e-8)));
	REQUIRE(approx(1e0) == linearWeight(flow, LatticePosition(0.5, 0.5, 1.0 - 1e-8)));
      }

        SECTION("testReadFromXML") {
            tinyxml2::XMLDocument doc;
            doc.Parse("<inlet>"
                      "  <normal units=\"dimensionless\" value=\"(0.0,1.0,1.0)\" />"
                      "  <position units=\"m\" value=\"(0.1,0.2,0.3)\" />"
                      "  <flowextension>"
                      "    <length units=\"m\" value=\"0.1\" />"
                      "    <radius units=\"m\" value=\"0.01\" />"
                      "    <fadelength units=\"m\" value=\"0.05\" />"
                      "  </flowextension>"
                      "</inlet>");
            auto converter = std::make_shared<util::UnitConverter>(0.5, 0.6, PhysicalPosition{0.7}, DEFAULT_FLUID_DENSITY_Kg_per_m3, 0.0);

            configuration::SimConfig conf;
            configuration::SimConfigReader reader("ignored.xml");
            configuration::CosinePressureIoletConfig iolet_conf;
            io::xml::Element inletEl = doc.FirstChildElement("inlet");
            reader.DoIOForBaseInOutlet(conf.GetSimInfo(), inletEl, iolet_conf);

            auto fe_conf = iolet_conf.flow_extension;
            REQUIRE(fe_conf.has_value());
            auto builder = UninitSimBuilder(conf, converter);
            auto flow = builder.BuildFlowExtension(*fe_conf);
	//auto flow = readFlowExtension(doc.FirstChildElement("inlet"), converter);
	// Normal is opposite direction compared to XML inlet definition
	REQUIRE(approx(-2e0 / std::sqrt(2)) == Dot(LatticePosition(0, 1, 1), flow->normal));
	auto const length = converter->ConvertToLatticeUnits("m", 0.1);
	REQUIRE(approx(length) == flow->length);
	REQUIRE(approx(converter->ConvertToLatticeUnits("m", 0.01)) == flow->radius);
	REQUIRE(approx(converter->ConvertToLatticeUnits("m", 0.05)) == flow->fadeLength);

	// position is at opposite end compared to XML inlet definition
	auto const position = converter->ConvertPositionToLatticeUnits(LatticePosition(0.1,
										      0.2,
										      0.3));
	REQUIRE(approx(position.x()) == flow->origin.x());
	REQUIRE(approx(position.y() + length / std::sqrt(2e0)) == flow->origin.y());
	REQUIRE(approx(position.z() + length / std::sqrt(2e0)) == flow->origin.z());

	// Check flow extension is positioned correctly
	auto isInFlow = [&flow, &converter](double x, double y, double z) {
	  LatticePosition const pos(x, y, z);
	  LatticePosition const conv = converter->ConvertPositionToLatticeUnits(pos);
	  return contains(flow, conv);
	};
	auto const l = 0.1 / std::sqrt(2);
	REQUIRE(isInFlow(0.1, 0.25, 0.35));
	REQUIRE(isInFlow(0.1, 0.2 + 0.001, 0.3 + 0.001));
	REQUIRE(isInFlow(0.1, 0.2 + l - 0.001, 0.3 + l - 0.001));
	REQUIRE(not isInFlow(0.1, 0.2 + l + 0.001, 0.3 + l + 0.001));
	REQUIRE(not isInFlow(0.1, 0.2 - 0.001, 0.3 - 0.001));
      }
    }

}

