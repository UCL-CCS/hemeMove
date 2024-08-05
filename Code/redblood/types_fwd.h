// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_TYPES_FWD_H
#define HEMELB_REDBLOOD_TYPES_FWD_H

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>

namespace hemelb::redblood
{
    class CellBase;

    namespace details
    {
      //! Stable comparison of cells across nodes
      struct CellUUIDComparison
      {
	inline bool operator()(std::shared_ptr<CellBase> const&a,
			std::shared_ptr<CellBase> const &b) const;
      };
    }
    //! Typical cell container type
    using CellContainer = std::set<std::shared_ptr<CellBase>, details::CellUUIDComparison>;
    //! \brief Container of template meshes
    //! \details An instance of this object is used to reference meshes across the simulation
    using TemplateCellContainer = std::map<std::string, std::shared_ptr<CellBase>> ;
    //! Function to insert cells somewhere
    using CellInserter = std::function<void(CellContainer::value_type)>;
    //! Type of callback for listening to changes to cells
    using CellChangeListener = std::function<void (const CellContainer &)>;
}
#endif

