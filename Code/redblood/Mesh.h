// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_MESH_H
#define HEMELB_REDBLOOD_MESH_H

#include <memory>
#include <array>
#include <vector>
#include <set>
#include <string>
#include <type_traits>
#include <map>

#include "redblood/MeshIdType.h"
#include "util/Vector3D.h"
#include "util/Matrix3D.h"
#include "util/UnitConverter.h"
#include "units.h"

class vtkPolyData;

namespace hemelb::redblood
{
    //! Holds raw mesh data
    //! Data is separated into vertices and triangular facets
    struct MeshData
    {
        //! Type of containers over indices
        using Facet = std::array<IdType, 3>;
        //! Facet container type
        using Facets = std::vector<Facet>;
        //! Vertex container type
        using Vertices = std::vector<LatticePosition>;
        //! Vertex container
        Vertices vertices;
        //! Facet container
        Facets facets;
    };
    static_assert(std::is_standard_layout_v<MeshData>, "Needed for MPI data type");
    static_assert(
        std::is_move_constructible_v<MeshData>
        and std::is_copy_constructible_v<MeshData>,
        "Explicit type characteristics"
    );

    LatticePosition barycentre(MeshData const &mesh);
    LatticePosition barycentre(MeshData::Vertices const &vertices);
    LatticeVolume volume(MeshData const &mesh);
    LatticeVolume volume(MeshData::Vertices const &vertices, MeshData::Facets const &facets);
    LatticeArea area(MeshData const &mesh);
    LatticeArea area(MeshData::Vertices const &vertices, MeshData::Facets const &facets);
    //! DEPRECATED. Orients facets outward, or inward. Algorithm cannot handle case of facet being coplanar with mesh barycentre.
    unsigned orientFacets(MeshData &mesh, bool outward = true);
    //! Orients facets inwards/outwards using VTK algorithm to determining outward facing direction. MeshData object should have been constructed from vtkPolyData object. See readMeshDataFromVTKPolyData.
    unsigned orientFacets(MeshData &mesh, vtkPolyData &polydata, bool outward = true);

    //! Holds raw topology data
    class MeshTopology
    {
      public:
        //! Type for map from vertices to facets
        using VertexToFacets = std::vector<std::set<IdType> >;
        //! Type for map from facets to its neighbors
        using FacetNeighbors = std::vector<std::array<IdType, 3> >;
        //! For each vertex, lists the facet indices
        VertexToFacets vertexToFacets;
        //! For each facet, lists the neighboring facets
        FacetNeighbors facetNeighbors;

        // Creates mesh topology from mesh data
        MeshTopology(MeshData const &mesh);
        // Empty topology
        MeshTopology()
        {
        }
    };
    static_assert(
        std::is_default_constructible_v<MeshTopology>
        and (not std::is_nothrow_default_constructible_v<MeshTopology>)
        and std::is_move_constructible_v<MeshTopology>
        and std::is_nothrow_move_constructible_v<MeshTopology>
        and std::is_copy_constructible_v<MeshTopology>
        and std::is_copy_assignable_v<MeshTopology>
        and (not std::is_nothrow_copy_assignable_v<MeshTopology>)
        and std::is_standard_layout_v<MeshTopology>
        and (not std::is_trivial_v<MeshTopology>),
        "Explicit type characteristics"
    );

    //! Triangular mesh
    class Mesh
    {
        //! Differentiates between two shallow and deep copy constructors
        struct deepcopy_tag
        {
        };

        //! Deep copy constructor
        Mesh(Mesh const &c, deepcopy_tag const &) :
            mesh(new MeshData(*c.mesh)), topology(new MeshTopology(*c.topology))
        {
        }

      public:
        //! Initializes mesh from mesh data
        Mesh(std::shared_ptr<MeshData> const &mesh) :
            mesh(mesh), topology(new MeshTopology(*mesh))
        {
        }
        //! Initializes mesh from mesh data and topology
        Mesh(std::shared_ptr<MeshData> const &mesh, std::shared_ptr<MeshTopology> const &topo) :
            mesh(mesh), topology(topo)
        {
        }
        //! Initialize mesh by copying data
        Mesh(MeshData const &data) :
            mesh(new MeshData(data)), topology(new MeshTopology(data))
        {
        }
        //! Shallow copy constructor
        Mesh(Mesh const &in) :
            mesh(in.mesh), topology(in.topology)
        {
        }

        //! Determines barycentre of mesh
        LatticePosition GetBarycentre() const
        {
          return barycentre(*mesh);
        }
        //! Computes volume of the mesh
        LatticeVolume GetVolume() const
        {
          return volume(*mesh);
        }
        //! Computes area of the mesh
        LatticeArea GetArea() const
        {
          return area(*mesh);
        }
        //! Connectivity data
        std::shared_ptr<const MeshTopology> GetTopology() const
        {
          return topology;
        }
        //! Returns mesh data
        std::shared_ptr<MeshData> GetData()
        {
          return mesh;
        }
        //! Returns mesh data
        std::shared_ptr<MeshData const> GetData() const
        {
          return mesh;
        }
        //! Returns mesh data
        void SetData(MeshData const &data)
        {
          mesh.reset(new MeshData(data));
        }

        //! Makes a copy of the mesh
        Mesh clone() const
        {
          return Mesh(*this, deepcopy_tag());
        }

        //! Scale mesh around barycentre
        void operator*=(Dimensionless const &scale);
        //! Scale by matrix around barycentre
        void operator*=(util::Matrix3D const &rotation);
        //! Translate mesh
        void operator+=(LatticePosition const &offset);
        //! Transform mesh
        void operator+=(std::vector<LatticePosition> const &displacements);

        //! Number of nodes
        size_t GetNumberOfNodes() const
        {
          return mesh->vertices.size();
        }
        //! Const access to vertices
        MeshData::Vertices &GetVertices()
        {
          return mesh->vertices;
        }
        //! Const access to vertices
        MeshData::Vertices const &GetVertices() const
        {
          return mesh->vertices;
        }
        //! Const access to vertices
        MeshData::Vertices::const_reference GetVertex(size_t site) const
        {
          return mesh->vertices[site];
        }
        //! Const access to facets
        MeshData::Facets const &GetFacets() const
        {
          return mesh->facets;
        }

        //! True if both mesh share the same data
        bool isSameData(Mesh const &other) const
        {
          return mesh == other.mesh and topology == other.topology;
        }

      protected:
        //! Holds actual data about the mesh
        std::shared_ptr<MeshData> mesh;
        //! Holds topology information
        std::shared_ptr<MeshTopology> topology;
    };
    static_assert(
        (not std::is_default_constructible_v<Mesh>)
        and (not std::is_nothrow_default_constructible_v<Mesh>)
        and std::is_move_constructible_v<Mesh>
        and (not std::is_nothrow_move_constructible_v<Mesh>)
        and std::is_copy_constructible_v<Mesh>
        and std::is_copy_assignable_v<Mesh>
        and std::is_nothrow_copy_assignable_v<Mesh>
        and std::is_standard_layout_v<Mesh>
        and (not std::is_trivial_v<Mesh>),
        "Explicit type characteristics"
    );

    //! Tetrahedron of a depth
    //! Depth refers to the number of triangular subdivision in each facet
    Mesh tetrahedron(unsigned int depth = 0);
    //! Flat triangular shape
    //! Depth refers to the number of triangular subdivision in each facet
    //! There are two initial facets, referring to the same three vertices.
    Mesh pancakeSamosa(unsigned int depth = 0);
    //! Refine a mesh by decomposing each facet into four triangles
    Mesh refine(Mesh mesh, unsigned int depth = 0);
    //! Creates an ico sphere (regular triangular mesh over a sphere)
    Mesh icoSphere(unsigned int depth = 0);

    //! Orients facet outward, or inward
    unsigned orientFacets(Mesh &mesh, bool outward = true);

}
#endif
