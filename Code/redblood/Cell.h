//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_PARTICLE_H
#define HEMELB_REDBLOOD_PARTICLE_H

#include "redblood/Mesh.h"
#include "redblood/Node2Node.h"
#include "units.h"

namespace hemelb { namespace redblood {

class CellBase {

  public:
  //! \brief Initializes mesh from mesh data
  //! \param [in] vertices: deformable vertices that define the cell. These
  //!    values are *not* modified by the scale.
  //! \param [in] origMesh: Original mesh. A shallow copy is made of this
  //!    object.
  //! \param [in] scale: scales template by a given amount
  //!    The scale is added during internal operations. The template will still
  //!    refer to the same data in memory.
  CellBase(
      MeshData::Vertices && vertices,
      Mesh const &origMesh,
      Dimensionless scale = 1e0
  ) : vertices_(std::move(vertices)), templateMesh(origMesh), scale_(scale) {
    assert(scale_ > 1e-12);
  }
  //! \brief Initializes mesh from mesh data
  //! \param [in] vertices: deformable vertices that define the cell. These
  //!    values are *not* modified by the scale.
  //! \param [in] origMesh: Original mesh. A shallow copy is made of this
  //!    object.
  //! \param [in] scale: scales template by a given amount
  //!    The scale is added during internal operations. The template will still
  //!    refer to the same data in memory.
  CellBase(
      MeshData::Vertices const & vertices,
      Mesh const &origMesh,
      Dimensionless scale = 1e0
  ) : vertices_(vertices), templateMesh(origMesh), scale_(scale) {
    assert(scale_ > 1e-12);
  }

  //! \brief Initializes mesh from mesh data
  //! \param [in] mesh: deformable vertices that define the cell are copied
  //!    from this mesh. These values are *not* modified by the scale.
  //! \param [in] origMesh: Original mesh. A shallow copy is made of this
  //!    object.
  //! \param [in] scale: scales template by a given amount
  //!    The scale is added during internal operations. The template will still
  //!    refer to the same data in memory.
  CellBase(Mesh const & mesh, Mesh const &origMesh, Dimensionless scale=1e0)
      : CellBase(mesh.GetVertices(), origMesh, scale) {}

  //! \brief Initializes mesh from mesh data
  //! \param [in] mesh: Modifyiable mesh and template. Deep copies are made of
  //!   both.
  CellBase(Mesh const & mesh)
      : CellBase(mesh.GetVertices(), mesh.clone()) {}

  //! \brief Initializes mesh from mesh data
  //! \param [in] mesh: Modifyiable mesh and template. Deep copies are made of
  //!   both
  CellBase(std::shared_ptr<MeshData> const & mesh)
       : CellBase(mesh->vertices, Mesh(*mesh)) {}

  void operator=(Mesh const& mesh) {
    templateMesh = mesh;
    vertices_ = templateMesh.GetVertices();
    scale_ = 1e0;
  }

  //! Because it is good practice
  virtual ~CellBase() {}

  //! Unmodified mesh
  Mesh const & GetTemplateMesh() const { return templateMesh; }
  //! Facets for the mesh
  MeshData::Facets const & GetFacets() const {
    return templateMesh.GetData()->facets;
  }
  //! Vertices of the cell
  MeshData::Vertices const & GetVertices() const { return vertices_; }
  //! Vertices of the cell
  MeshData::Vertices& GetVertices() { return vertices_; }
  //! Topology of the (template) mesh
  std::shared_ptr<MeshTopology const> GetTopology() const {
    return templateMesh.GetTopology();
  }
  size_t GetNumberOfNodes() const { return vertices_.size(); }

  //! Facet bending energy
  virtual PhysicalEnergy operator()() const = 0;
  //! Facet bending energy
  virtual PhysicalEnergy operator()(
      std::vector<LatticeForceVector> &in) const = 0;

  //! Scale mesh around barycenter
  void operator*=(Dimensionless const &scale);
  //! Scale by matrix around barycenter
  void operator*=(util::Matrix3D const &scale);
  //! Translate mesh
  void operator+=(LatticePosition const &offset);
  //! Transform mesh
  void operator+=(std::vector<LatticePosition> const &displacements);

  MeshData::Vertices::value_type GetBarycenter() const {
    return barycenter(vertices_);
  }

  //! Scale to apply to the template mesh
  void SetScale(Dimensionless scale) {
    assert(scale > 1e-12);
    scale_ = scale;
  }
  //! Scale to apply to the template mesh
  Dimensionless GetScale() const { return scale_; }

  protected:
   //! Holds list of vertices for this cell
   MeshData::Vertices vertices_;
   //! Unmodified original mesh
   Mesh templateMesh;
   //! Scale factor for the template;
   Dimensionless scale_;
};

//! Deformable cell for which energy and forces can be computed
class Cell : public CellBase {
public:
  //! Holds all physical parameters
  struct Moduli {
    //! Bending energy parameter
    PhysicalPressure bending;
    //! Surface energy parameter
    PhysicalPressure surface;
    //! Surface volume parameter
    PhysicalPressure volume;
    //! Skalak dilation modulus
    PhysicalPressure dilation;
    //! Skalak strain modulus
    PhysicalPressure strain;
    Moduli() : bending(0), surface(0), volume(0), dilation(0), strain(0) {};
  } moduli;
  //! Node-wall interaction
  Node2NodeForce nodeWall;

  // inheriting constructors
  using CellBase::CellBase;

  //! Facet bending energy
  virtual PhysicalEnergy operator()() const override;
  //! Facet bending energy
  virtual PhysicalEnergy operator()(
      std::vector<LatticeForceVector> &in) const override;

private:
  // Computes facet bending energy over all facets
  PhysicalEnergy facetBending_() const;
  // Computes facet bending energy over all facets
  PhysicalEnergy facetBending_(std::vector<LatticeForceVector> &forces) const;
};

//! Typical cell container type
typedef std::vector<std::shared_ptr<CellBase>> CellContainer;

}} // namespace hemelb::redblood
#endif
