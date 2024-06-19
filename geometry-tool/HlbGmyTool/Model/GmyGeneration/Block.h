// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_BLOCK_H
#define HLBGMYTOOL_GMY_BLOCK_H

#include <vtkSmartPointer.h>
#include <vector>

class vtkOBBTree;

#include "Domain.h"
#include "Index.h"
namespace hemelb::gmytool::gmy {

class Site;
using SiteVec = std::vector<Site*>;
using SiteIterator = SiteVec::iterator;

class Block {
 public:
  // typedef std::vector<Site*>::iterator iterator;

  Block(Domain&, const Index&, const unsigned int&);
  ~Block();

  Site& GetGlobalSite(const Index&);
  Site& GetLocalSite(const Index&);

  inline SiteIterator begin() { return this->sites.begin(); }

  inline SiteIterator end() { return this->sites.end(); }

  inline Domain& GetDomain() const { return this->domain; }
  inline const Index& GetIndex() const { return this->index; }
  vtkSmartPointer<vtkOBBTree> CreateOBBTreeModel(double extraSize) const;

  const Site& Middle() const { return *sites[sites.size() / 2]; }

 protected:
  unsigned int size;
  const Index index;
  const Index min;
  const Index max;
  Domain& domain;
  SiteVec sites;

  inline unsigned int TranslateIndex(const Index& ind) {
    return (ind[0] * this->size + ind[1]) * this->size + ind[2];
  }
  friend class NeighbourIteratorBase;
  friend class LaterNeighbourIterator;
};

}  // namespace hemelb::gmytool::gmy
#endif  // HLBGMYTOOL_GMY_BLOCK_H
