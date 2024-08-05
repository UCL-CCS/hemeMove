// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLBASE_H
#define HEMELB_REDBLOOD_CELLBASE_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <set>
#include <utility>

#include "redblood/Mesh.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
#   ifndef NDEBUG
    //! \brief Function checks characteristics of implementation class CellBase::CellData
    //! \data Nothing more than a bunch of static asserts. CellData is a private member, so we
    //! need a friend somewhere to check its characteristics.
    void checkCellDataCharacteristics();
#   endif
    class CellBase
    {
#     ifndef NDEBUG
        friend void checkCellDataCharacteristics();
#     endif
      public:
        //! \brief Initializes mesh from mesh data
        //! \param [in] verticesIn: deformable vertices that define the cell. These
        //!    values are *not* modified by the scale.
        //! \param [in] origMesh: Original mesh. A shallow copy is made of this
        //!    object.
        //! \param [in] scale: scales template by a given amount
        //!    The scale is added during internal operations. The template will still
        //!    refer to the same data in memory.
        CellBase(MeshData::Vertices &&verticesIn, Mesh const &origMesh, LatticeDistance scale = 1e0,
                 std::string const &templateName = "default");
        //! \brief Initializes mesh from mesh data
        //! \param [in] verticesIn: deformable vertices that define the cell. These
        //!    values are *not* modified by the scale.
        //! \param [in] origMesh: Original mesh. A shallow copy is made of this
        //!    object.
        //! \param [in] scale: scales template by a given amount
        //!    The scale is added during internal operations. The template will still
        //!    refer to the same data in memory.
        CellBase(MeshData::Vertices const &verticesIn, Mesh const &origMesh, LatticeDistance scale =
                     1e0,
                 std::string const &templateName = "default");

        //! \brief Initializes mesh from mesh data
        //! \param [in] mesh: deformable vertices that define the cell are copied
        //!    from this mesh. These values are *not* modified by the scale.
        //! \param [in] origMesh: Original mesh. A shallow copy is made of this
        //!    object.
        //! \param [in] scale: scales template by a given amount
        //!    The scale is added during internal operations. The template will still
        //!    refer to the same data in memory.
        CellBase(Mesh const &mesh, Mesh const &origMesh, LatticeDistance scale = 1e0,
                 std::string const &templateName = "default");

        //! \brief Initializes mesh from mesh data
        //! \param [in] mesh: Modifyiable mesh and template. Deep copies are made of
        //!   both.
        CellBase(Mesh const &mesh);
        //! \brief Initializes mesh from mesh data
        //! \param [in] mesh: Modifyiable mesh and template. Deep copies are made of
        //!   both
        CellBase(std::shared_ptr<MeshData> const &mesh);
        //! Copy constructor
        //! References same template mesh
        CellBase(CellBase const &cell);

        //! Tag to choose shallow clone constructor
        struct shallow_clone
        {
        };
        //! Shallow copy constructor
        //! References same data
        CellBase(CellBase const &cell, shallow_clone const&) :
            data(cell.data)
        {
        }
        //! Because it is good practice
        virtual ~CellBase()
        {
        }

        void operator=(Mesh const &mesh);

        //! Unmodified mesh
        Mesh const &GetTemplateMesh() const;
        //! Unmodified mesh
        std::string const &GetTemplateName() const;
        //! Modifies template cell
        void SetTemplateName(std::string const&);
        //! Facets for the mesh
        MeshData::Facets const &GetFacets() const;
        //! Vertices of the cell
        MeshData::Vertices const &GetVertices() const;
        //! Vertices of the cell
        MeshData::Vertices &GetVertices();
        //! Topology of the (template) mesh
        std::shared_ptr<MeshTopology const> GetTopology() const;
        site_t GetNumberOfNodes() const;

        //! Facet bending energy
        virtual LatticeEnergy operator()() const = 0;
        //! Facet bending energy - pretty printing for shared ptr
        LatticeEnergy Energy() const
        {
          return operator()();
        }
        //! Facet bending energy
        virtual LatticeEnergy operator()(std::vector<LatticeForceVector> &in) const = 0;
        //! Facet bending energy - pretty printing for shared ptr
        LatticeEnergy Energy(std::vector<LatticeForceVector> &in) const
        {
          return operator()(in);
        }

        //! Scale mesh around barycentre
        void operator*=(Dimensionless const &);
        //! Linear transform of each vertex, centered around barycentre
        void operator*=(util::Matrix3D const &);
        //! Translate mesh
        void operator+=(LatticePosition const &offset);
        //! Translate mesh
        void operator-=(LatticePosition const &offset)
        {
          return operator+=(-offset);
        }
        //! Transform mesh
        void operator+=(std::vector<LatticePosition> const &displacements);

        MeshData::Vertices::value_type GetBarycentre() const;
        LatticeVolume GetVolume() const
        {
          return volume(GetVertices(), GetTemplateMesh().GetFacets());
        }

        //! Scale to apply to the template mesh
        void SetScale(LatticeDistance scale);
        //! Scale to apply to the template mesh
        LatticeDistance GetScale() const;

        // cloneImpl is virtual and returns a pointer to abstract class
        // clone will be overriden. It will call cloneImpl and cast it to derived type.
        // The hoop jumping is necessary since we are returning managed pointers.
        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<CellBase> clone() const
        {
          return cloneImpl();
        }

        //! A unique identifier for the cell
        boost::uuids::uuid const & GetTag() const;
        //! Sets the cell's unique identifier
        void SetTag(boost::uuids::uuid const & uuid);

        //! Computes average edge length of cell
        double GetAverageEdgeLength() const;

      protected:
        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<CellBase> virtual cloneImpl() const = 0;
        //! allows separation of data and behaviors
        class CellData;
        //! Holds data
        std::shared_ptr<CellData> data;

        //! Allows more control by derived classes
        CellBase(std::shared_ptr<CellData> const &dataIn) :
            data(dataIn)
        {
        }
    };
    static_assert(
        (not std::is_default_constructible_v<CellBase>)
        and (not std::is_nothrow_default_constructible_v<CellBase>)
        and (not std::is_move_constructible_v<CellBase>)
        and (not std::is_nothrow_move_constructible_v<CellBase>)
        and (not std::is_copy_constructible_v<CellBase>)
        and std::is_copy_assignable_v<CellBase>
        and std::is_nothrow_copy_assignable_v<CellBase>
        and (not std::is_standard_layout_v<CellBase>)
        and (not std::is_trivial_v<CellBase>),
        "Explicit type characteristics"
    );

    // Holds the data of CellBase
    // Although its own data members are public, it can only be accessed by CellBase and derived
    // classes. So encapsulation is not broken. Indeed, the definition of this class should not be
    // available to anyone but the CellBase class family. Adding Get/Setters would mean more code
    // for the exact same functionality.
    class CellBase::CellData
    {
      public:
        CellData(MeshData::Vertices &&verticesIn, Mesh const &origMesh, LatticeDistance scaleIn =
                     1e0,
                 std::string templateName = "default") :
            vertices(std::move(verticesIn)), templateMesh(origMesh), scale(scaleIn),
                tag(boost::uuids::random_generator()()), templateName(templateName)
        {
          assert(scale > 1e-12);
        }
        CellData(MeshData::Vertices const &verticesIn, Mesh const &origMesh,
                 LatticeDistance scaleIn = 1e0, std::string const &templateName = "default") :
            vertices(verticesIn), templateMesh(origMesh), scale(scaleIn),
                tag(boost::uuids::random_generator()()), templateName(templateName)
        {
          assert(scale > 1e-12);
        }
        CellData(CellData const& c) :
            vertices(c.vertices), templateMesh(c.templateMesh), scale(c.scale),
                tag(boost::uuids::random_generator()()), templateName(c.templateName)
        {
        }
        CellData(CellData && c) :
            vertices(std::move(c.vertices)), templateMesh(std::move(c.templateMesh)),
                scale(c.scale), tag(std::move(c.tag)), templateName(std::move(c.templateName))
        {
        }
        // Constructor that can copy the tag
        CellData(Mesh const &origMesh, boost::uuids::uuid const &uuid,
                 std::string const &templateName) :
            templateMesh(origMesh), scale(1e0), tag(uuid), templateName(templateName)
        {
          assert(scale > 1e-12);
        }
        //! Holds list of vertices for this cell
        MeshData::Vertices vertices;
        //! Unmodified original mesh
        Mesh templateMesh;
        //! Scale factor for the template;
        LatticeDistance scale;
        //! Uuid tag
        boost::uuids::uuid tag;
        //! \brief tag of the cell from which this one was cloned
        //! \details In practice, all cells are generated from a few templates. This tag identifies
        //! that template. If not given on input, then it is set to "default"
        std::string templateName;

    };
  }
} // namespace hemelb::redblood
#endif
