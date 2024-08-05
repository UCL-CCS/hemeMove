// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <fstream>
#include <numeric>
#include "redblood/CellBase.h"

namespace hemelb
{
  namespace redblood
  {
    CellBase::CellBase(MeshData::Vertices &&verticesIn, Mesh const &origMesh, Dimensionless scaleIn,
                       std::string const & templateName) :
        data(new CellData(std::move(verticesIn), origMesh, scaleIn, templateName))
    {
    }
    CellBase::CellBase(MeshData::Vertices const &verticesIn, Mesh const &origMesh,
                       Dimensionless scaleIn, std::string const & templateName) :
        data(new CellData(verticesIn, origMesh, scaleIn, templateName))
    {
    }
    CellBase::CellBase(CellBase const& cell) :
        data(new CellData(*cell.data))
    {
    }
    CellBase::CellBase(Mesh const &mesh, Mesh const &origMesh, Dimensionless scaleIn,
                       std::string const & templateName) :
        data(new CellData(mesh.GetVertices(), origMesh, scaleIn, templateName))
    {
    }
    CellBase::CellBase(Mesh const &mesh) :
        data(new CellData(mesh.GetVertices(), mesh.clone()))
    {
    }
    CellBase::CellBase(std::shared_ptr<MeshData> const &mesh) :
        data(new CellData(mesh->vertices, Mesh(*mesh)))
    {
    }

    void CellBase::operator=(Mesh const &mesh)
    {
      data->templateMesh = mesh;
      data->vertices = data->templateMesh.GetVertices();
      data->scale = 1e0;
    }

    //! Unmodified mesh
    Mesh const &CellBase::GetTemplateMesh() const
    {
      return data->templateMesh;
    }
    std::string const & CellBase::GetTemplateName() const
    {
      return data->templateName;
    }
    void CellBase::SetTemplateName(std::string const& name)
    {
      data->templateName = name;
    }
    //! Facets for the mesh
    MeshData::Facets const &CellBase::GetFacets() const
    {
      return data->templateMesh.GetData()->facets;
    }
    //! Vertices of the cell
    MeshData::Vertices const &CellBase::GetVertices() const
    {
      return data->vertices;
    }
    //! Vertices of the cell
    MeshData::Vertices &CellBase::GetVertices()
    {
      return data->vertices;
    }
    //! Topology of the (template) mesh
    std::shared_ptr<MeshTopology const> CellBase::GetTopology() const
    {
      return data->templateMesh.GetTopology();
    }
    site_t CellBase::GetNumberOfNodes() const
    {
      return static_cast<site_t>(data->vertices.size());
    }
    MeshData::Vertices::value_type CellBase::GetBarycentre() const
    {
      return barycentre(data->vertices);
    }

    //! Scale to apply to the template mesh
    void CellBase::SetScale(Dimensionless scaleIn)
    {
      assert(data->scale > 1e-12);
      data->scale = scaleIn;
    }
    //! Scale to apply to the template mesh
    Dimensionless CellBase::GetScale() const
    {
      return data->scale;
    }

    void CellBase::operator*=(Dimensionless const &scaleIn)
    {
      auto const barycentre = GetBarycentre();

      for (auto &vertex : data->vertices)
      {
        vertex = (vertex - barycentre) * scaleIn + barycentre;
      }
    }
    void CellBase::operator*=(util::Matrix3D const &rotation)
    {
      auto const barycentre = GetBarycentre();

      for (auto &vertex : data->vertices)
      {
        rotation.timesVector(vertex - barycentre, vertex);
        vertex += barycentre;
      }
    }
    void CellBase::operator+=(LatticePosition const &offset)
    {
      for (auto &vertex : data->vertices)
      {
        vertex += offset;
      }
    }
    void CellBase::operator+=(std::vector<LatticePosition> const &displacements)
    {
      assert(displacements.size() == data->vertices.size());
      auto i_disp = displacements.begin();

      for (auto &vertex : data->vertices)
      {
        vertex += * (i_disp++);
      }
    }

    boost::uuids::uuid const & CellBase::GetTag() const
    {
      return data->tag;
    }
    void CellBase::SetTag(boost::uuids::uuid const & uuid)
    {
      data->tag = uuid;
    }

    double CellBase::GetAverageEdgeLength() const
    {
      std::vector<double> edgeLengths;

      for (auto &facet : data->templateMesh.GetFacets())
      {
        edgeLengths.push_back( (data->vertices[facet[0]] - data->vertices[facet[1]]).GetMagnitude());
        edgeLengths.push_back( (data->vertices[facet[1]] - data->vertices[facet[2]]).GetMagnitude());
        edgeLengths.push_back( (data->vertices[facet[2]] - data->vertices[facet[0]]).GetMagnitude());
      }

      return std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0) / edgeLengths.size();
    }

#   ifndef NDEBUG
    void checkCellDataCharacteristics()
    {
      static_assert(
          (not std::is_default_constructible<CellBase::CellData>::value)
          and (not std::is_nothrow_default_constructible<CellBase::CellData>::value)
          and std::is_move_constructible<CellBase::CellData>::value
          and (not std::is_nothrow_move_constructible<CellBase::CellData>::value)
          and std::is_copy_constructible<CellBase::CellData>::value
          and (not std::is_copy_assignable<CellBase::CellData>::value)
          and (not std::is_nothrow_copy_assignable<CellBase::CellData>::value)
          and std::is_standard_layout<CellBase::CellData>::value
          and (not std::is_pod<CellBase::CellData>::value),
          "Explicit type characteristics"
      );
    }
#   endif
  }
} // hemelb::redblood

