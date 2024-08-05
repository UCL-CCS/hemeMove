// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "io/xml.h"
#include "configuration/CommandLine.h"
#include "configuration/SimBuilder.h"
#include "configuration/SimConfigReader.h"
#include "redblood/CellControllerBuilder.h"
#include "redblood/MeshIO.h"
#include "redblood/CellIO.h"
#include "redblood/FaderCell.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/SimConfBuildHelp.h"

template <typename To, typename From>
auto dynamic_unique_cast(std::unique_ptr<From>&& src) {
    From* src_raw = src.release();
    To* dest_raw = dynamic_cast<To*>(src_raw);
    if (dest_raw) {
        // Worked!
    } else {
        // Fail - need to do something with the resource, guess we assign back to src
        src.reset(src_raw);
    }
    return std::unique_ptr<To>(dest_raw);
}

namespace hemelb::tests
{
    using namespace redblood;
    namespace fs = std::filesystem;

    TEST_CASE_METHOD(helpers::FolderTestFixture, "CellIOTests", "[redblood]") {
        //tinyxml2::XMLDocument doc;
        io::xml::Document doc;
        LatticeDistance scale;

        auto converter = std::make_shared<util::UnitConverter>(0.5, 0.6, LatticePosition(1, 2, 3), 1000.0, 0.0);

        CopyResourceToTempdir("red_blood_cell.txt");
        CopyResourceToTempdir("empty_for_relative_paths.xml");
        configuration::SimConfig config;
        configuration::SimConfigReader reader("empty_for_relative_paths.xml");
        //auto config = std::make_unique<UninitSimConfig>("empty_for_relative_paths.xml");
        auto builder = UninitSimBuilder(config, converter);
        auto cell_builder = CellControllerBuilder(converter);
        {
            // This allocates an element, but doesn't add to the document.
            // (Note: elements are owned by their document)
            auto parent = doc.AddChild("parent");
            auto cell = parent.AddChild("cell");

            auto shape = cell.AddChild("shape");
            shape.SetAttribute("mesh_path", "red_blood_cell.txt");
            shape.SetAttribute("mesh_format", "Krueger");

            auto scaleXML = cell.AddChild("scale");
            scale = 1.5;
            scaleXML.SetAttribute("value", scale);
            scaleXML.SetAttribute("units", "m");

            auto moduli = cell.AddChild("moduli");
            auto add_stuff = [&moduli](char const *name, char const *units, Dimensionless value) {
                auto elem = moduli.AddChild(name);
                elem.SetAttribute("value", value);
                elem.SetAttribute("units", units);
            };
            add_stuff("surface", "lattice", 2e0);
            add_stuff("dilation", "lattice", 0.58);
            add_stuff("bending", "Nm", 2e-18);
        }

        auto approx = Approx(0.0).margin(1e-12);

        // Reads cell with minimum stuff
        SECTION("testReadCellWithDefaults") {
            // Remove moduli, so we get default behavior
            auto cellEl = doc.GetRoot().GetChildOrThrow("cell");
            cellEl.GetChildOrThrow("moduli").Delete();

            auto cellConf = reader.readCell(cellEl);
            auto cell = dynamic_unique_cast < Cell const>(cell_builder.build_cell(cellConf));
            REQUIRE(cell);
            auto const kruegerIO = redblood::KruegerMeshIO{};
            auto const data = kruegerIO.readFile("red_blood_cell.txt", true);
            REQUIRE(static_cast<site_t>(data->vertices.size()) == cell->GetNumberOfNodes());
            REQUIRE(approx(converter->ConvertToLatticeUnits("m", scale)) == cell->GetScale());
            REQUIRE(approx(1e0) == cell->moduli.volume);
            REQUIRE(approx(1e0) == cell->moduli.surface);
            REQUIRE(approx(converter->ConvertToLatticeUnits("Nm", 2e-19)) == cell->moduli.bending);
            REQUIRE(approx(0.75) == cell->moduli.dilation);
            REQUIRE(approx(converter->ConvertToLatticeUnits("N/m", 5e-6)) == cell->moduli.strain);
        }

        SECTION("testReadCellModuli") {
            auto cellEl = doc.GetRoot().GetChildOrThrow("cell");
            auto cellConf = reader.readCell(cellEl);
            auto const moduli = cell_builder.build_cell_moduli(cellConf.moduli);
            REQUIRE(approx(1e0) == moduli.volume);
            REQUIRE(approx(2e0) == moduli.surface);
            REQUIRE(approx(converter->ConvertToLatticeUnits("Nm", 2e-18)) == moduli.bending);
            REQUIRE(approx(0.58) == moduli.dilation);
            REQUIRE(approx(converter->ConvertToLatticeUnits("N/m", 5e-6)) == moduli.strain);
        }

        SECTION("testReadMeshTemplates") {
            const char *xml_text = "<parent>"
                                   "  <inlets>"
                                   "   <inlet>"
                                   "     <condition type=\"pressure\" subtype=\"cosine\">"
                                   "       <amplitude value=\"0\" units=\"Pa\" />"
                                   "       <mean value=\"0\" units=\"Pa\" />"
                                   "       <phase value=\"0\" units=\"rad\" />"
                                   "       <period value=\"1\" units=\"s\" />"
                                   "     </condition>"
                                   "     <normal units=\"dimensionless\" value=\"(0.0,0.0,1.0)\" />"
                                   "     <position units=\"m\" value=\"(0.0,0.0,-0.024)\" />"
                                   "     <flowextension>"
                                   "       <length units=\"m\" value=\"0.1\" />"
                                   "       <radius units=\"m\" value=\"0.01\" />"
                                   "       <fadelength units=\"m\" value=\"0.05\" />"
                                   "     </flowextension>"
                                   "   </inlet>"
                                   "  </inlets>"
                                   "  <outlets>"
                                   "    <outlet>"
                                   "     <condition type=\"pressure\" subtype=\"cosine\">"
                                   "       <amplitude value=\"0\" units=\"Pa\" />"
                                   "       <mean value=\"0\" units=\"Pa\" />"
                                   "       <phase value=\"0\" units=\"rad\" />"
                                   "       <period value=\"1\" units=\"s\" />"
                                   "     </condition>"
                                   "      <normal units=\"dimensionless\" value=\"(0.0,0.0,-1.0)\" />"
                                   "      <position units=\"m\" value=\"(0.0,0.0,0.024)\" />"
                                   "      <flowextension>"
                                   "        <length units=\"m\" value=\"0.1\" />"
                                   "        <radius units=\"m\" value=\"0.01\" />"
                                   "        <fadelength units=\"m\" value=\"0.05\" />"
                                   "      </flowextension>"
                                   "    </outlet>"
                                   "  </outlets>"
                                   "  <redbloodcells>"
                                   "    <templates>"
                                   "      <cell>"
                                   "        <shape mesh_path=\"red_blood_cell.txt\" mesh_format=\"Krueger\" />"
                                   "        <scale units=\"m\" value=\"0.6\"/>"
                                   "      </cell>"
                                   "     <cell name=\"joe\">"
                                   "       <shape mesh_path=\"red_blood_cell.txt\" mesh_format=\"Krueger\" />"
                                   "       <scale units=\"m\" value=\"0.5\"/>"
                                   "     </cell>"
                                   "   </templates>"
                                   "  </redbloodcells>"
                                   "</parent>";
            io::xml::Document document;
            document.LoadString(xml_text);
            auto root = document.GetRoot();
            auto cellsEl = root.GetChildOrThrow("redbloodcells").GetChildOrThrow("templates");
            auto tc_conf = reader.readTemplateCells(cellsEl);
            configuration::GlobalSimInfo sim_info;
            auto in_conf = reader.DoIOForInOutlets(sim_info, root.GetChildOrThrow("inlets"));
            auto out_conf = reader.DoIOForInOutlets(sim_info, root.GetChildOrThrow("outlets"));
            auto ins = builder.BuildIolets(in_conf);
            auto outs = builder.BuildIolets(out_conf);

            auto cells = cell_builder.build_template_cells(
                    tc_conf, configuration::MakeCountedIoletView(ins), configuration::MakeCountedIoletView(outs)
            );
            REQUIRE(size_t(2) == cells->size());
            REQUIRE(size_t(1) == cells->count("default"));
            REQUIRE(size_t(1) == cells->count("joe"));
            auto const default_ = std::static_pointer_cast<FaderCell>((*cells)["default"]);
            auto const joe = std::static_pointer_cast<FaderCell>((*cells)["joe"]);
            REQUIRE(default_->GetTemplateName() == "default");
            REQUIRE(joe->GetTemplateName() == "joe");
            REQUIRE(Approx(converter->ConvertToLatticeUnits("m", 0.6)).margin(1e-8)
                    == default_->GetScale());
            REQUIRE(Approx(converter->ConvertToLatticeUnits("m", 0.5)).margin(1e-8)
                    == joe->GetScale());
            REQUIRE(static_cast<bool>(default_->GetIOlets()));
            REQUIRE(static_cast<bool>(joe->GetIOlets()));
            REQUIRE(default_->GetIOlets() == joe->GetIOlets());
            REQUIRE(size_t(2) == joe->GetIOlets()->size());
        }

        // Reads cell with minimum stuff
        SECTION("ReadWriteReadCell") {
            // Remove moduli, so we get default behavior
            auto cellEl = doc.GetRoot().GetChildOrThrow("cell");
            cellEl.GetChildOrThrow("moduli").Delete();

            auto cellConf = reader.readCell(cellEl);
            auto cell = std::shared_ptr(cell_builder.build_cell(cellConf));
            {
                // Write the cell mesh and barycentre data
                auto cells = CellContainer{};
                cells.insert(cell);
                auto outConf = configuration::CellOutputConfig{10, false};
                configuration::CommandLine cmdline({"hemelb", "-in", "empty_for_relative_paths.xml"});
                auto pathmgr = std::make_shared<configuration::PathManager>(cmdline, Comms().OnIORank(), Comms().Size());
                auto simState = std::make_shared<lb::SimulationState>(0.1, 1000);
                auto full_out = cell_builder.build_full_cell_output(outConf, simState, pathmgr, Comms());
                auto summ_out = cell_builder.build_summary_cell_output(outConf, simState, pathmgr, Comms());

                full_out(cells);
                summ_out(cells);
            }
            {
                // Check barycentre
                auto bary_path = fs::path{"results/Cells/0/barycentres.rbc"};
                REQUIRE(fs::exists(bary_path));
                auto bci = CellBarycentreInput{bary_path};
                auto ncells = bci.ReadHeader();
                REQUIRE(ncells == 1);
                auto cell_summary = bci.ReadRows(Comms(), 0, 1);
                auto [uuid, bcpos] = cell_summary[0];
                REQUIRE(uuid == cell->GetTag());
                REQUIRE(bcpos == ApproxV(cell->GetBarycentre()));

                // Check mesh
                char tag[36+5];
                to_chars(uuid, tag);
                std::strncpy(tag + 36, ".vtp", 5);
                auto mesh_path = fs::path{"results/Cells/0"} / tag;
                REQUIRE(fs::exists(mesh_path));

                auto mesh_reader = VTKMeshIO{};
                auto actual_mesh = mesh_reader.readFile(mesh_path, false);
                REQUIRE(cell->GetVertices() == actual_mesh->vertices);
                REQUIRE(cell->GetFacets() == actual_mesh->facets);
            }
        }
    }
}
