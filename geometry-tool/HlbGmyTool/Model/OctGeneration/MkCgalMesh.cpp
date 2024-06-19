// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "MkCgalMesh.h"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace hemelb::gmytool::oct {

class SurfaceCreator : public CGAL::Modifier_base<CgalPolyhedron::HalfedgeDS> {
 public:
  typedef CgalPolyhedron::HalfedgeDS HDS;

  SurfaceCreator(const std::vector<Vector>& ptsIn,
                 const std::vector<Index>& polysIn)
      : points(ptsIn), triangles(polysIn) {}
  void operator()(HDS& hds);

 private:
  const std::vector<Vector>& points;
  const std::vector<Index>& triangles;
};

void SurfaceCreator::operator()(HDS& hds) {
  CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
  B.begin_surface(points.size(), triangles.size(), 0);

  typedef typename HDS::Vertex Vertex;
  typedef typename Vertex::Point Point;

  for (auto pt : points)
    B.add_vertex(Point(pt.x, pt.y, pt.z));

  for (size_t i = 0; i < triangles.size(); ++i) {
    auto& tri = triangles[i];
    // VTK polygons can contain lines where two vertexes are identical. Forget
    // these
    if ((tri[0] != tri[1]) & (tri[0] != tri[2]) & (tri[1] != tri[2])) {
      auto face = B.begin_facet();
      B.add_vertex_to_facet(tri[0]);
      B.add_vertex_to_facet(tri[1]);
      B.add_vertex_to_facet(tri[2]);
      B.end_facet();
      // Set the id to the original cell's id.
      face->id() = i;
    } else {
      std::cout << "Eliminated degenerate vertex: " << tri << std::endl;
    }
  }

  B.end_surface();
  B.remove_unconnected_vertices();
}

CgalMeshPtr MkCgalMesh(const std::vector<Vector>& points,
                       const std::vector<Index>& triangles) {
  SurfaceCreator modifier(points, triangles);
  auto mesh = std::make_shared<CgalPolyhedron>();
  mesh->delegate(modifier);
  return mesh;
}

}  // namespace hemelb::gmytool::oct
