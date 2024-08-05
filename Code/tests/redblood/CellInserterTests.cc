// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "configuration/SimConfigReader.h"
#include "io/xml.h"
#include "redblood/RBCInserter.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/SimConfBuildHelp.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb::tests
{
    using namespace io::xml;
    using namespace redblood;

    TEST_CASE_METHOD(FlowExtensionFixture, "CellInserterTests", "[redblood]") {
        auto converter = std::make_shared<util::UnitConverter>(0.5, 0.6, PhysicalPosition::Zero(), 1000.0, 0.0);
        LatticeTime every = 10, offset = 5;
        TemplateCellContainer cells;
        cells.emplace("joe", std::make_shared<Cell>(tetrahedron()));

        auto getDocument = [&](LatticeDistance radius = 1e0, int numInserters = 0) {
            std::ostringstream sstr;
            sstr << "<hemelbsettings>\n<inlets>\n<inlet>\n"
                    "  <condition type=\"pressure\" subtype=\"file\">\n"
                    "  <path value=\"ignored.dat\" />\n"
                    "  </condition>\n"
                    "  <normal units=\"dimensionless\" value=\"(0.0,1.0,1.0)\" />\n"
                    "  <position units=\"m\" value=\"(0.1,0.2,0.3)\" />\n"
                    "  <flowextension>\n"
                    "    <length units=\"m\" value=\"0.5\" />\n"
                    "    <radius units=\"m\" value=\"" << radius << "\" />\n"
                    "    <fadelength units=\"m\" value=\"0.4\" />\n"
                    "  </flowextension>\n";
            for (int i = 0; i < numInserters; ++i)
            {
                sstr << "  <insertcell template=\"joe\">\n"
                     << "    <every units=\"s\" value=\"" << every << "\"/>\n"
                     << "    <offset units=\"s\" value=\"" << offset << "\"/>\n"
                     << "    <seed value=\"854" << i << "\" />\n"
                     << "  </insertcell>\n";
            }
            sstr << "</inlet>\n</inlets>\n</hemelbsettings>\n";
            Document doc;
            doc.LoadString(sstr.str().c_str());
            return doc;
        };

        auto readRBCInserters = [&](io::xml::Element const& inletsNode, TemplateCellContainer const& templateCells) {
            using namespace configuration;
            SimConfig conf;
            SimConfigReader reader("ignored.xml");
            auto inlet_confs = reader.DoIOForInOutlets(conf.GetSimInfo(), inletsNode);
            auto builder = UninitSimBuilder(conf, converter);
            auto inlets = builder.BuildIolets(inlet_confs);

            auto cell_builder = redblood::CellControllerBuilder(converter);
            return cell_builder.build_cell_inserters(inlet_confs, configuration::MakeCountedIoletView(inlets), templateCells);
        };

        SECTION("testNoPeriodicInsertion") {
            auto doc = getDocument(1e0, 0);
            *cells["joe"] *= 0.1e0;
            auto const inserter = readRBCInserters(doc.GetRoot().GetChildOrThrow("inlets"),
                                                   cells);
            REQUIRE(not inserter);
        }

        SECTION("testCellOutsideFlowExtension") {
            auto doc = getDocument(1e0, 1);
            REQUIRE_THROWS_AS(readRBCInserters(doc.GetRoot().GetChildOrThrow("inlets"),
                                               cells),
                              hemelb::Exception);
        }

        SECTION("testPeriodicInsertion") {
            // Creates an inserter and checks it exists
            auto doc = getDocument(1, 2);
            *cells["joe"] *= 0.1e0;
            auto const inserter = readRBCInserters(doc.GetRoot().GetChildOrThrow("inlets"),
                                                   cells);
            REQUIRE(inserter);

	// all calls up to offset result in node added cell
	int num_calls = 0;
	auto addCell = [&num_calls](CellContainer::value_type const&)
	  {
	    // This function is called by inserter with a new cell
	    // It is meant to actually do the job of adding cells to the simulation
	    ++num_calls;
	  };
	auto const dt = converter->ConvertTimeToPhysicalUnits(1e0);
	auto const N = static_cast<int>(std::floor(offset / dt));
	for (int i(0); i < N; ++i) {
	  inserter(addCell);
	  REQUIRE(0 == num_calls);
	}

	// Now first call should result in adding a cell
	inserter(addCell);
	REQUIRE(2 == num_calls);
	num_calls = 0;

	// Now should not output a cell until "every" time has gone by
	// Do it thrice for good measure
	for (int k = 0; k < 3; ++k) {
	  for (int i(0); i < int(std::floor(every / dt)) - 1; ++i) {
	    inserter(addCell);
	    REQUIRE(0 == num_calls);
	  }
	  // And should call it!
	  inserter(addCell);
	  REQUIRE(2 == num_calls);
	  num_calls = 0;
	}
        }

      SECTION("testTranslation") {
	auto const identity = rotationMatrix(LatticePosition(0, 0, 1),
					     LatticePosition(0, 0, 1));
	RBCInserterWithPerturbation inserter([]()
					     { return true;},
					     cells["joe"]->clone(), identity, 0e0, 0e0, LatticePosition(2,
													0,
													0),
					     LatticePosition(0, 4, 0));

	auto const barycentre = cells["joe"]->GetBarycentre();
	for (size_t i(0); i < 500; ++i) {
	  auto const cell = inserter.drop();
	  auto const n = cell->GetBarycentre();
	  REQUIRE(std::abs(n.x() - barycentre.x()) <= 2e0);
	  REQUIRE(std::abs(n.y() - barycentre.y()) <= 4e0);
	  REQUIRE(Approx(barycentre.z()).margin(1e-8) == n.z());
	}
      }

      SECTION("testSpatialOffset") {
	auto doc = getDocument(1e0, 1);
	helpers::ModifyXMLInput(doc,
				{ "inlets", "inlet", "insertcell", "offset", "value" },
				0.0);
	helpers::ModifyXMLInput(doc,
				{ "inlets", "inlet", "insertcell", "every", "value" },
				0.4);
	*cells["joe"] *= 0.1e0;
	CellContainer::value_type current_cell;
	auto addCell = [&current_cell](CellContainer::value_type const& cell) {
	  current_cell = cell;
	};
    // Have to copy as doc is consumed by reading
    auto doc1 = doc.DeepCopy();
	auto const inserter = readRBCInserters(doc1.GetRoot().GetChildOrThrow("inlets"),
                                           cells);
	REQUIRE(inserter);
	inserter(addCell);
	REQUIRE(current_cell);
	auto const cell0 = current_cell;

	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "x", "units" }, "m");
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "x", "value" }, 0.1);
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "y", "units" }, "m");
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "y", "value" }, 0.1);
    auto doc2 = doc.DeepCopy();
	auto const insertTranslated = readRBCInserters(doc2.GetRoot().GetChildOrThrow("inlets"),
                                                   cells);
	REQUIRE(insertTranslated);
	insertTranslated(addCell);
	REQUIRE(current_cell);
	auto const trans = cell0->GetBarycentre() - current_cell->GetBarycentre();
	auto approx = Approx(0.0).margin(1e-8);
	REQUIRE(approx(0) == Dot(trans, util::Vector3D{0, 1, 1}));
	REQUIRE(approx(0.1 / 0.6 * 0.1 / 0.6 * 2e0) == trans.GetMagnitudeSquared());

	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "z", "units" }, "m");
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "z", "value" }, 0.1);
	auto const insertWithZ = readRBCInserters(doc.GetRoot().GetChildOrThrow("inlets"),
                                              cells);
	REQUIRE(insertWithZ);
	insertWithZ(addCell);
	REQUIRE(current_cell);
	auto const transZ = cell0->GetBarycentre() - current_cell->GetBarycentre();
	REQUIRE(approx(0.1 / 0.6 * std::sqrt(2)) == Dot(transZ, util::Vector3D{0, 1, 1}));
	REQUIRE(approx(0.1 / 0.6 * 0.1 / 0.6 * 3e0) == transZ.GetMagnitudeSquared());
      }

    }

  }
