// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_BUILDCGALPOLYGON_H
#define HLBGMYTOOL_GMY_BUILDCGALPOLYGON_H

#include "CGALtypedef.h"

#include "Block.h"

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

class vtkPoints;
class vtkCellArray;
class vtkIntArray;

namespace hemelb::gmytool::gmy {

template <class HDS>
class BuildCGALPolygon : public CGAL::Modifier_base<HDS> {
 public:
  BuildCGALPolygon(vtkPoints* ptsin,
                   vtkCellArray* polysin,
                   vtkIntArray* IoletIdArrayIn) {
    this->pts = ptsin;
    this->polys = polysin;
    this->IoletIdArray = IoletIdArrayIn;
  }
  void operator()(HDS& hds);

 private:
  vtkPoints* pts;
  vtkCellArray* polys;
  vtkIntArray* IoletIdArray;
};

}  // namespace hemelb::gmytool::gmy

#endif  // HLBGMYTOOL_GMY_BUILDCGALPOLYGON_H
