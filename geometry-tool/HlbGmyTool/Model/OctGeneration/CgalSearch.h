// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_CGALSEARCH_H
#define HLBGMYTOOL_OCT_CGALSEARCH_H

#include "Cgal.h"

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
//#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

namespace hemelb::gmytool::oct {

// typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel, CgalPolyhedron>
// Primitive;
typedef CGAL::AABB_face_graph_triangle_primitive<CgalPolyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABBTraits;

typedef CGAL::AABB_tree<AABBTraits> CgalSearchTree;

}  // namespace hemelb::gmytool::oct
#endif
