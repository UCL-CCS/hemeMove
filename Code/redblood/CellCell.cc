//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <vector>
#include "redblood/CellCell.h"
#include "redblood/Cell.h"
#include "redblood/Interpolation.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace redblood
  {
    namespace
    {
      template<class T>
      size_t figureNearness(DivideConquer<T> &dnc, LatticeVector const &key,
                         LatticePosition const &vertex, LatticeDistance const &haloLength)
      {
        if (haloLength + haloLength > dnc.GetBoxSize())
        {
          return static_cast<size_t>(Borders::CENTER);
        }

        size_t result = 0;

        for (size_t d(1); d < (1 << 6); d <<= 1)
        {
          LatticePosition const translated(direction<LatticePosition::value_type>(d) * haloLength);

          if (not (key == dnc.DowngradeKey(vertex + translated)))
          {
            result |= d;
          }
        }

        return result bitor static_cast<size_t>(Borders::CENTER);
      }

      template<class T>
      CellReference initCellRef(DivideConquer<T> &dnc, CellContainer::const_iterator cellid,
                                site_t nodeid, LatticeVector const &key,
                                LatticePosition const &vertex, LatticeDistance const &haloLength)
      {
        return {cellid, nodeid, figureNearness(dnc, key, vertex, haloLength)};
      }

      void initializeCells(DivideConquer<CellReference> &dnc, MeshData::Vertices const &vertices,
                           CellContainer::const_iterator cellid, LatticeDistance haloLength)
      {
        typedef DivideConquer<CellReference> DnC;
        typedef DnC::key_type key_type;
        typedef MeshData::Vertices::const_iterator vertex_iterator;
        vertex_iterator i_first = vertices.begin();
        vertex_iterator const i_end = vertices.end();

        for (site_t i(0); i_first != i_end; ++i_first, ++i)
        {
          key_type const key = dnc.DowngradeKey(*i_first);
          dnc.insert(key, initCellRef(dnc, cellid, i, key, *i_first, haloLength));
        }
      }

      // avoids a warning
#     ifndef HEMELB_DOING_UNITTESTS
      void initializeCells(DivideConquer<CellReference> &dnc, CellContainer const &cells,
                           LatticeDistance haloLength)
      {
        CellContainer::const_iterator i_first = cells.begin();
        CellContainer::const_iterator const i_end = cells.end();

        for (site_t i(0); i_first != i_end; ++i_first, ++i)
        {
          initializeCells(dnc, (*i_first)->GetVertices(), i_first, haloLength);
        }
      }
#     endif

      // Compare distance between vertices
      template<class T_FUNCTION>
      bool nextDistance(T_FUNCTION const &cellOrdering, DivideConquerCells::const_iterator &first,
                        DivideConquerCells::const_iterator const &end,
                        DivideConquerCells::const_iterator const &main, LatticeDistance dist)
      {
        auto const mainCell = main.GetCell();
        typedef DivideConquerCells::const_iterator cit;
        auto goodCellPair = [&mainCell, &cellOrdering](cit const &i)
        {
          return cellOrdering(i.GetCell(), mainCell);
        };
        auto goodDistance = [&main, &dist](cit const &i)
        {
          return (*main - *i).GetMagnitude() < dist;
        };
        for (; first != end; ++first)
        {
          if (goodCellPair(first) and goodDistance(first))
          {
            return true;
          }
        }
        return false;
      }
    } // anonymous namespace
#ifndef HEMELB_DOING_UNITTESTS
    //! Constructor
    DivideConquerCells::DivideConquerCells(CellContainer const &cells, LatticeDistance boxsize,
                                           LatticeDistance halosize) :
        DivideConquer<CellReference>(boxsize), haloLength(halosize), cells(cells)
    {
      initializeCells(*static_cast<base_type *>(this), GetCells(), haloLength);
    }

    void DivideConquerCells::update()
    {
      iterator i_first = begin();
      iterator const i_end = end();

      for (; i_first != i_end;)
      {
        key_type const key = base_type::DowngradeKey(*i_first);
        i_first.GetCellReference().nearBorder = figureNearness(*this, key, *i_first, haloLength);

        if (not (key == i_first.GetKey()))
        {
          base_type::insert(key, i_first.GetCellReference());
          base_type::iterator const next_iter = base_type::erase((base_type::iterator) i_first);
          i_first = iterator(*this, next_iter);
        }
        else
        {
          ++i_first;
        }
      }
    }

    DivideConquerCells::const_range DivideConquerCells::operator()(LatticeVector const &pos) const
    {
      base_type::const_range const boxrange = base_type::equal_range(pos);
      return const_range(const_iterator(*this, boxrange.first),
                         const_iterator(*this, boxrange.second));
    }

    bool DivideConquerCells::pair_range::nextDist()
    {
      typedef decltype(owner.cells) Cells;
      typedef Cells::const_reference Input;
      auto strictly_less = Cells::key_compare();
      auto cellOrdering = [&strictly_less](Input _a, Input _b)
      {
        return _a != _b and not strictly_less(_a, _b);
      };
      return nextDistance(cellOrdering, currents.second, ends.second, currents.first, maxdist);
    }

    bool DivideConquerCells::pair_range::doBox()
    {
      LatticeVector const key(currents.first.GetKey() + *box_iterator);
      DivideConquerCells::const_range const boxits = owner(key);

      if (box_iterator->x == 0 and box_iterator->y == 0 and box_iterator->z ==0)
      {
        currents.second = currents.first;
        ++currents.second;
      }
      else
      {
        currents.second = boxits.first;
      }

      ends.second = boxits.second;
      return nextDist();
    }

    bool DivideConquerCells::pair_range::operator++()
    {
      if (not is_valid())
      {
        return false;
      }

      // First try and finds next pair in current range
      if (currents.second != ends.second)
      {
        ++currents.second;

        if (nextDist())
        {
          return true;
        }
      }

      // If reaches here, then go to next box
      for(++box_iterator; box_iterator; ++box_iterator)
      {
        if(doBox())
        {
          return true;
        }
      }

      // If reaches here, then should increment main iterator and start with same box
      if (++currents.first == ends.first)
      {
        return false;
      }

      box_iterator = BorderBoxIterator(currents.first.GetNearBorder());
      return doBox() ? true : operator++();
    }

    DivideConquerCells::pair_range::pair_range(DivideConquerCells const &owner,
                                               iterator const &begin, iterator const &end,
                                               LatticeDistance maxdist) :
        maxdist(maxdist), box_iterator(begin.GetNearBorder()),
        currents(begin, end), ends(end, end), owner(owner)
    {
      // No throw garantee. Makes iterator invalid instead.
      try
      {
        // Could be invalid from start
        if (not is_valid())
        {
          return;
        }

        // Iterates to first valid item, if any
        if (not doBox())
        {
          operator++();
        }
      }
      catch (std::exception const &e)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("*** Encountered error while initializing pair iterator: %s\n",
                                                      e.what());
        currents.first = ends.first;
      }
      catch (...)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("*** Encountered error while initializing pair iterator.");
        currents.first = ends.first;
      }
    }

    DivideConquerCells::pair_range DivideConquerCells::pair_begin(LatticeDistance maxdist) const
    {
      return pair_range(*this, begin(), end(), maxdist);
    }

    bool DivideConquerCells::insert(CellContainer::value_type cell)
    {
      auto inserted = cells.insert(cell);
      if (inserted.second)
        initializeCells(*this, (*inserted.first)->GetVertices(), inserted.first, haloLength);
      return true;
    }

    bool DivideConquerCells::remove(CellContainer::value_type cell)
    {
      auto const cellIterator = cells.find(cell);
      if (cellIterator == cells.end())
        return false;
      auto i_first = base_type::cbegin();
      auto const i_end = base_type::cend();
      while (i_first != i_end)
      {
        auto const i_current = i_first;
        // Move to next element *before* possibly erasing current element
        // Not sure whether i_first can be incremented after an erase.
        ++i_first;
        if (*i_current->second.cellIterator == cell)
        {
          base_type::erase(i_current);
        }
      }
      cells.erase(cellIterator);
      return true;
    }
#endif
  }
} // namespace hemelb::redblood
