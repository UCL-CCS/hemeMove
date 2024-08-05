// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <sstream>
#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "resources/Resource.h"
#include "util/UnitConverter.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb::tests
{
    using namespace redblood;

    TEST_CASE("RedBloodMeshDataIOTests") {
        redblood::KruegerMeshIO io = {};
        std::string filename = resources::Resource("red_blood_cell.txt").Path();
        std::shared_ptr<MeshData> mesh = io.readFile(filename, true);
      
        SECTION("testReadMesh") {
            REQUIRE(mesh);
            REQUIRE(mesh->vertices.size() == 12);
            REQUIRE(mesh->facets.size() == 20);

            using Vector3d = util::Vector3D<double> ;
            Vector3d vfirst(0.850650808352040, 0.525731112119134, 0.000000000000000);
            Vector3d vlast(0, -0.767597772408033, -0.319039022201661);
            REQUIRE(mesh->vertices.front() == ApproxV(vfirst));
            REQUIRE(mesh->vertices.back() == ApproxV(vlast));

            // Is any of the vertices of the facet equal to value?
            auto any = [](MeshData::Facet const &vec, MeshData::Facet::value_type value) -> bool {
                return vec[0] == value or vec[1] == value or vec[2] == value;
            };

            REQUIRE(any(mesh->facets.front(), 0));
            REQUIRE(any(mesh->facets.front(), 8));
            REQUIRE(any(mesh->facets.front(), 4));
            REQUIRE(any(mesh->facets.back(), 11));
            REQUIRE(any(mesh->facets.back(), 7));
            REQUIRE(any(mesh->facets.back(), 5));
        }

        SECTION("testWriteMesh") {
            //auto conv = util::UnitConverter(1, 1, LatticePosition(0, 0, 0), 1000.0, 0.0);
            auto output = io.writeString(*mesh, nullptr);
            std::shared_ptr<MeshData> other = io.readString(output, true);
            REQUIRE(other->vertices.size() == mesh->vertices.size());
            REQUIRE(other->facets.size() == mesh->facets.size());
            REQUIRE(ApproxV(mesh->vertices.front()) == other->vertices.front());
            REQUIRE(ApproxV(mesh->vertices.back()) == other->vertices.back());

            for (size_t i = 0; i < 3; ++i)
            {
                REQUIRE(mesh->facets.front()[i] == other->facets.front()[i]);
                REQUIRE(mesh->facets.back()[i] == other->facets.back()[i]);
            }
        }

    }
}
