// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iterator>
#include "redblood/buffer/Columns.h"
#include "util/Matrix3D.h"

namespace hemelb::redblood::buffer
{
    namespace detail
    {
        LatticeDistance maxExtension(MeshData::Vertices const &vertices,
                                     LatticePosition const &direction)
        {
          auto const barycentre = redblood::barycentre(vertices);
          LatticeDistance result(0);
          auto const normalised = direction.GetNormalised();
          auto maxdist = [&normalised, &result, &barycentre](LatticePosition const &b)
          {
            result = std::max(result, std::abs(Dot(normalised, b - barycentre)));
          };
          std::for_each(vertices.begin(), vertices.end(), maxdist);
          return 2e0 * result;
        }

        LatticePosition maxExtensions(MeshData::Vertices const &vertices,
                                      LatticePosition const &col, LatticePosition const& normal)
        {
          LatticeDistance const x = maxExtension(vertices, normal);
          LatticeDistance const y = maxExtension(vertices, Cross(col, normal));
          LatticeDistance const z = maxExtension(vertices, col);
          return {x, y, z};
        }
      }

    ColumnPositionIterator::ColumnPositionIterator(std::shared_ptr<Cylinder const> i_cylinder,
                                                   MeshData::Vertices const& vertices,
                                                   LatticePosition cellAxis,
                                                   LatticePosition colAxis,
                                                   LatticeDistance separation) :
            IndexIterator(LatticeVector(0, 0, 0), LatticeVector(0, 0, 0)), cylinder(std::move(i_cylinder))
    {
        if (cellAxis.GetMagnitude() < 1e-8)
        {
          throw Exception() << "Cell axis cannot be 0";
        }
        if (colAxis.GetMagnitude() < 1e-8)
        {
          throw Exception() << "Column axis cannot be 0";
        }
        if (std::abs(Dot(colAxis, cylinder->normal)) > 1e-8)
        {
          throw Exception() << "Column axis should be perpendicular to cylinder axis";
        }

        // Rotation is opposite to the one that will be applied to the mesh
        auto const antiRot = rotationMatrix(colAxis, cellAxis);
        auto const extents = detail::maxExtensions(vertices, antiRot * colAxis, antiRot * cylinder->normal)
	  + LatticePosition{separation};

        max = LatticeVector(
                std::numeric_limits<LatticeCoordinate>::max(),
                std::ceil(cylinder->radius / extents.y()),
                std::ceil(cylinder->radius / extents.z())
        );
        min = LatticeVector(0, -max.y(), -max.z());
        current = min;

        major = colAxis.GetNormalised() * extents.z();
        minor = Cross(colAxis.GetNormalised(), cylinder->normal.GetNormalised()) * extents.y();
        depth = cylinder->normal.GetNormalised() * extents.x();
        spacing = (extents + LatticePosition{separation}) / 2.0;

        if (not IsInside())
          operator++();
      }

      bool ColumnPositionIterator::IsInside() const
      {
        // avoids missing defs for operator^
        LatticeDistance const a = major.GetMagnitudeSquared()
            / std::pow(cylinder->radius - spacing.z(), 2e0);
        LatticeDistance const b = minor.GetMagnitudeSquared()
            / std::pow(cylinder->radius - spacing.y(), 2e0);
        return (current.z() * current.z()) * a + (current.y() * current.y()) * b < 1e0;
      }

      void ColumnPositionIterator::operator++()
      {
        assert(IsValid());
        IndexIterator::operator++();

#       ifndef NDEBUG
        size_t breaker = 100, iter = 0;
#       endif
        while (not IsInside())
        {
          IndexIterator::operator++();
#         ifndef NDEBUG
          if ( (++iter) >= breaker)
          {
            throw Exception() << "Probably in an infinite loop";
          }
#         endif
        }
      }

      ColumnCellDrop::ColumnCellDrop(std::shared_ptr<Cylinder const> cylinder,
                                     CellContainer::value_type cell, LatticePosition cellAxis,
                                     LatticePosition colAxis, LatticeDistance wallSeparation) :
          iterator(cylinder, cell->GetVertices(), cellAxis, colAxis, wallSeparation),
              templateCell(cell->clone())
      {
        // Cell is rotated to correct orientation
        *templateCell *= rotationMatrix(cellAxis, colAxis);
        // And centered at zero
        *templateCell -= templateCell->GetBarycentre();
      }

      CellContainer::value_type ColumnCellDrop::operator()()
      {
        ++iterator;
        std::shared_ptr<CellBase> result(templateCell->clone());
        *result += *iterator;
        return result;
      }

} // hemelb::redblood
