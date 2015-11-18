//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#include "redblood/VertexBag.h"
#include "redblood/Interpolation.h"

namespace hemelb
{
  namespace redblood
  {
    namespace
    {
      // avoids a warning
#     ifndef HEMELB_DOING_UNITTESTS
      proc_t get_proc(geometry::LatticeData const &latDat, LatticePosition const &pos)
      {
        proc_t procid;
        site_t siteid;
        return latDat.GetContiguousSiteId(pos, procid, siteid) ?
          procid :
          std::numeric_limits<proc_t>::max();
      }
#     endif
      //! Set of procs affected by this position
      //! \param[in] latDat will tell us which site belongs to which proc
      //! \param[in] iterator a  stencil iterator going over affected lattice points
      template<class STENCIL, class T_FUNC> std::set<proc_t> procsAffectedByPosition(
          T_FUNC get_proc, InterpolationIterator<STENCIL> &&iterator, proc_t avoid)
      {
        std::set<proc_t> result;
        for (; iterator.IsValid(); ++iterator)
        {
          auto const procid = get_proc(*iterator);
          if (procid != avoid)
          {
            result.insert(procid);
          }
        }
        return result;
      }
      //! Set of procs affected by this position
      //! \param[in] latDat will tell us which site belongs to which proc
      //! \param[in] position for which to figure out affected processes
      //! \param[in] stencil giving interaction range
      template<class STENCIL, class T_FUNC> std::set<proc_t> procsAffectedByPosition(
          T_FUNC get_proc, LatticePosition const &position,
          proc_t avoid = std::numeric_limits<proc_t>::max())
      {
        return procsAffectedByPosition(get_proc, interpolationIterator<STENCIL>(position), avoid);
      }

      template<class STENCIL, class T_FUNC>
      std::map<size_t, std::shared_ptr<VertexBag>> splitVertices(
          T_FUNC get_proc, std::shared_ptr<CellBase const> cell,
          proc_t avoid = std::numeric_limits<proc_t>::max())
      {
        std::map<size_t, std::shared_ptr<VertexBag>> result;
        for (auto const &vertex : cell->GetVertices())
        {
          auto const regions = procsAffectedByPosition<STENCIL>(get_proc, vertex, avoid);
          for (auto const region : regions)
          {
            auto const i_bag = result.find(region);
            if (i_bag == result.end())
            {
              result.emplace(region, std::make_shared<VertexBag>(cell, vertex));
            }
            else
            {
              i_bag->second->addVertex(vertex);
            }
          }
        }
        return result;
      }

    }

#   ifndef HEMELB_DOING_UNITTESTS
    VertexBag::VertexBag(std::shared_ptr<CellBase const> parent) :
            CellBase(std::shared_ptr<CellBase::CellData>(new CellData(parent->GetTemplateMesh(),
                                                                      parent->GetTag(),
                                                                      parent->GetTemplateName())))
    {
    }
    VertexBag::VertexBag(std::shared_ptr<CellBase const> parent, LatticePosition vertex) :
            CellBase(std::shared_ptr<CellBase::CellData>(new CellData(parent->GetTemplateMesh(),
                                                                      parent->GetTag(),
                                                                      parent->GetTemplateName())))
    {
      addVertex(vertex);
    }

    VertexBag::VertexBag(boost::uuids::uuid const &tag, std::string const &templateName)
      : CellBase(
          std::make_shared<CellBase::CellData>(
            Mesh(std::make_shared<MeshData>(MeshData())), tag, templateName))
    {
    }

    template <class STENCIL>
    std::map<size_t, std::shared_ptr<VertexBag>> splitVertices(
        std::shared_ptr<CellBase const> cell, geometry::LatticeData const &latticeData,
        proc_t selfRegion)
    {
      auto proc_getter = [&latticeData](LatticePosition const &position)
      {
        return get_proc(latticeData, position);
      };
      return splitVertices<STENCIL>(proc_getter, cell, selfRegion);
    }
#   endif
  }
} // hemelb::redblood
