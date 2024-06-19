// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "PolyDataGenerator.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdList.h>
#include <vtkMatrix4x4.h>
#include <vtkOBBTree.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkXMLPolyDataWriter.h>

#include "InconsistentFluidnessError.h"
#include "Neighbours.h"
#include "Site.h"

#include "BuildCGALPolygon.cpp"
#include "BuildCGALPolygon.h"
#include "CGALtypedef.h"
// because the class is templated we need to include the cpp file to avoid link
// errors. the alternative is a explicit instantiation in BuildCGALPolygon.cpp
// i.e add template class BuildCGALPolygon<HalfedgeDS>; however that seems to be
// slower.

#include "Debug.h"

#include "io/formats/geometry.h"

#include "Block.h"

namespace hemelb::gmytool::gmy {

using namespace io::formats;

PolyDataGenerator::PolyDataGenerator()
    : GeometryGenerator(), ClippedSurface(NULL) {
  this->Locator = vtkOBBTree::New();
  // this->Locator->SetNumberOfCellsPerNode(32); // the default
  this->Locator->SetTolerance(1e-9);
  this->hitPoints = vtkPoints::New();
  this->hitCellIds = vtkIdList::New();
}

PolyDataGenerator::~PolyDataGenerator() {}

void PolyDataGenerator::ComputeBounds(double bounds[6]) const {
  this->ClippedSurface->GetBounds(bounds);
}

void PolyDataGenerator::ClosePolygon(void) {
  Halfedge_iterator j;
  Halfedge_iterator k;
  Halfedge_iterator l;
  Halfedge_iterator m;
  Halfedge_handle newedge;

  cout << "The polyhedron has " << ClippedCGALSurface->size_of_facets()
       << " facets " << ClippedCGALSurface->size_of_halfedges() << " halfedges "
       << ClippedCGALSurface->size_of_border_halfedges() << " border halfedges "
       << ClippedCGALSurface->size_of_vertices() << " vertices " << endl;
  int ncompremoved = ClippedCGALSurface->keep_largest_connected_components(1);
  cout << "Removing " << ncompremoved << " non connected components" << endl;
  cout << "The polyhedron has " << ClippedCGALSurface->size_of_facets()
       << " facets " << ClippedCGALSurface->size_of_halfedges() << " halfedges "
       << ClippedCGALSurface->size_of_border_halfedges() << " border halfedges "
       << ClippedCGALSurface->size_of_vertices() << " vertices " << endl;

  while (!this->ClippedCGALSurface->is_closed()) {
    for (j = ClippedCGALSurface->border_halfedges_begin();
         j != ClippedCGALSurface->halfedges_end(); ++j) {
      // We itterate over all border half edges. Half of them are opposite
      // border edegs and thus not real borders so test this first. it.next() it
      // the next border edge along the hole. This might not be the same as ++it
      k = j->next();
      l = k->next();
      m = l->next();
      if (j->is_border() && l->is_border() && k->is_border()) {
        // First find size of hole and make sure that this is simple. i.e.
        // No two different halfedges points to the same vertex.
        // If this is the case we delete one of them and start over.
        Halfedge_iterator tempit = k;
        Halfedge_iterator tempit2;
        int sizehole = 1;
        bool Simplehole = true;
        while (tempit != j) {
          for (tempit2 = j; tempit2 != tempit; tempit2 = tempit2->next()) {
            if (tempit2->vertex() == tempit->vertex()) {
              Simplehole = false;
              ClippedCGALSurface->erase_facet(tempit2->opposite());
              break;  // remove one facet and start over.
            }
          }
          if (!Simplehole) {
            break;  // ugly Break inermost while loop.
          }
          tempit = tempit->next();
          ++sizehole;
        }

        if (Simplehole) {
          // cout << "Hole has " << sizehole << endl;
          if (j != m) {  // more than 3 edges in hole. Have to subdivide
            newedge = ClippedCGALSurface->add_facet_to_border(j, l);
            // cout << "Filling " << j->vertex()->point() <<  " , " <<
            // k->vertex()->point() << " to " << newedge->vertex()->point() <<
            // endl;
          } else {
            // cout << "Closing " << j->vertex()->point() << " to " <<
            // k->vertex()->point() << endl;
            newedge = ClippedCGALSurface->fill_hole(j);
          }
          // now find the iolet type from neighbors.
          tempit = newedge;
          bool start = true;
          std::vector<size_t> neighborIds;
          while (tempit != newedge || start) {
            start = false;
            if (!tempit->opposite()->is_border()) {
              neighborIds.push_back(tempit->opposite()->facet()->id());
            }
            tempit = tempit->next();
          }
          std::vector<size_t>::iterator max =
              std::max_element(neighborIds.begin(), neighborIds.end());
          std::vector<size_t>::iterator min =
              std::min_element(neighborIds.begin(), neighborIds.end());
          if (*max != *min) {  // if max is min they are all the same.
            throw GenerationErrorMessage(
                "Could not determin IOlet id of inserted facet.");
          } else {
            newedge->facet()->id() = *max;
          }
          break;  // we have to normalize the border and restart the forloop
                  // since
          // filling with add_facet_to_border invalidated the order of
          // border/non border edges.
        } else {
          // cout << "Could not fill since hole is not simple, deleted connected
          // facet instead" << endl;
          break;
        }
      }
    }
    ClippedCGALSurface->normalize_border();
  }
}

void PolyDataGenerator::CreateCGALPolygon(void) {
  vtkPoints* pts;
  vtkCellArray* polys;
  std::clock_t start;
  double duration;
  start = std::clock();

  polys = this->ClippedSurface->GetPolys();
  pts = this->ClippedSurface->GetPoints();
  vtkIntArray* Iolets = this->IoletIdArray;
  this->triangle =
      std::make_unique<BuildCGALPolygon<HalfedgeDS>>(pts, polys, Iolets);
  this->ClippedCGALSurface = std::make_unique<Polyhedron>();
  this->ClippedCGALSurface->delegate(*this->triangle);
  this->ClippedCGALSurface->normalize_border();
  if (!this->ClippedCGALSurface->is_closed()) {
    cout << "The polygon is not closed most likely due to non manifold edges "
            "ignored. Will atempt to close it."
         << endl;
    ClosePolygon();
  }

  if (!this->ClippedCGALSurface->is_closed()) {
    throw GenerationErrorMessage("Created surface is not closed.");
  } else if (!ClippedCGALSurface->is_pure_triangle()) {
    throw GenerationErrorMessage(
        "Created surface is not pure triangles, cannot voxelize.");
  } else if (!ClippedCGALSurface->is_valid()) {
    throw GenerationErrorMessage(
        "Created surface is not valid, cannot voxelize.");
  } else {
    cout << "Succesfully created closed polygon from input" << endl;
    cout << "The polyhedron has " << ClippedCGALSurface->size_of_facets()
         << " facets " << ClippedCGALSurface->size_of_halfedges()
         << " halfedges " << ClippedCGALSurface->size_of_border_halfedges()
         << " border halfedges " << ClippedCGALSurface->size_of_vertices()
         << " vertices " << endl;
    bool write_out = true;
    if (write_out) {
      std::ofstream out;
      out.open("exportedsurface.off");
      out << (*ClippedCGALSurface);
      out.close();
    }
  }
  this->AABBtree =
      std::make_unique<Tree>(this->ClippedCGALSurface->facets_begin(),
                             this->ClippedCGALSurface->facets_end());
  duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
  std::cout << "Preprocessing took: " << duration << " s " << endl;
}

void PolyDataGenerator::PreExecute(void) {
  // Build our locator.
  this->Locator->SetDataSet(this->ClippedSurface);
  this->Locator->BuildLocator();

  /*
   * Get the scalars associated with the surface polygons. The scalars hold
   * the index of the Iolet which they represent, -1 meaning they aren't an
   * Iolet.
   */
  this->IoletIdArray = vtkIntArray::SafeDownCast(
      this->ClippedSurface->GetCellData()->GetScalars());
  this->CreateCGALPolygon();
  if (this->IoletIdArray == NULL) {
    throw GenerationErrorMessage(
        "Error getting Iolet ID array from clipped surface");
  }
}

/*
 * Given a site with known fluidness, examine the links to not-yet-visited
 * neighbouring sites. If the neighbours have unknown fluidness, set that.
 *
 * Since we wish to examine each link only once, this will set link properties
 * from neighbour => site as well as site => neighbour
 *
 */
void PolyDataGenerator::ClassifySite(Site& site) {
  if (!site.IsFluidKnown) {
    throw GenerationErrorMessage("The start site is not known cannot continue");
  }
  for (LaterNeighbourIterator neighIt = site.begin(); neighIt != site.end();
       ++neighIt) {
    Site& neigh = *neighIt;
    unsigned int iNeigh = neighIt.GetNeighbourIndex();
    int nHits = Intersect(site, neigh);
    // Four cases: fluid-fluid, solid-solid, fluid-solid and solid-fluid.
    // Will handle the last two together.
    if (site.IsFluid == neigh.IsFluid) {
      if (site.IsFluid) {
        // Fluid-fluid, must set CutType::NONE for both
        site.Links[iNeigh].Type = geometry::CutType::NONE;
        neigh.Links[Neighbours::inverses[iNeigh]].Type =
            geometry::CutType::NONE;
      } else {
        // solid-solid, nothing to do.
      }
    } else {
      // They differ, figure out which is fluid and which is solid.
      Site* fluid;
      Site* solid;

      // Index of the solid site from the fluid site.
      int iSolid;
      // Index of the point in this->hitPoints we're considering as the
      // hit of interest (i.e. the one closest to the fluid site).
      int iHit;

      Vector hitPoint;
      PointCGAL hitPointCGAL;
      SegmentCGAL hitsegmentCGAL;
      int hitCellId;
      double distancetol = 0.01;
      std::vector<int> CloseIoletIDS;
      if (nHits == -1) {  // unclassified intersection need to be carefull.
        if (site.IsFluid) {
          fluid = &site;
          solid = &neigh;
          iSolid = iNeigh;
          int n = 0;
          iHit = n;
          // here we iterate over all intersections until we
          // find a wall or the rest are futher away than the tol.
          for (std::vector<Object_Primitive_and_distance>::iterator distit =
                   IntersectionCGAL.begin();
               distit != IntersectionCGAL.end(); ++distit) {
            if (distit->second > IntersectionCGAL[0].second + distancetol) {
              iHit = n;
              break;
            }
            int tempioletId = distit->first.second->id() - 2;
            // shifting back from unsigned
            if (tempioletId < 0) {
              break;  // hit a wall no need to continue.
            }  // what if we hit both an inlet and outlet. Can that ever happen?
            ++n;
          }

        } else {
          fluid = &neigh;
          solid = &site;
          iSolid = Neighbours::inverses[iNeigh];
          int n = IntersectionCGAL.size() - 1;
          iHit = n;
          for (std::vector<Object_Primitive_and_distance>::reverse_iterator
                   distit = IntersectionCGAL.rbegin();
               distit != IntersectionCGAL.rend(); ++distit) {
            if (distit->second < IntersectionCGAL.back().second - distancetol) {
              iHit = n;
              break;  // ignoring the following intersections, they are to far
                      // away.
            }
            int tempioletId =
                distit->first.second->id() - 2;  // shifting back from unsigned
            if (tempioletId < 0) {
              break;  // hit a wall no need to continue.
            }  // what if we hit both an inlet and outlet. Can that ever happen?
            --n;
          }
        }
      } else {  // normal intersection just take the closest point.

        if (site.IsFluid) {
          fluid = &site;
          solid = &neigh;
          iSolid = iNeigh;
          iHit = 0;
        } else {
          fluid = &neigh;
          solid = &site;
          iSolid = Neighbours::inverses[iNeigh];
          iHit = nHits - 1;
        }
      }
      Object_Primitive_and_distance hitpoint_triangle_dist =
          IntersectionCGAL[iHit];

      if (CGAL::assign(hitPointCGAL, hitpoint_triangle_dist.first.first)) {
        // we do an explicite cast to double here.
        // The cast to double is only needed if we use an exact_construction
        // kernel. Otherwise this is already a double but keeping this in makes
        // it posible to change the kernel for testing.
        hitPoint = Vector(CGAL::to_double(hitPointCGAL.x()),
                          CGAL::to_double(hitPointCGAL.y()),
                          CGAL::to_double(hitPointCGAL.z()));
      } else if (CGAL::assign(hitsegmentCGAL,
                              hitpoint_triangle_dist.first.first)) {
        hitPointCGAL =
            CGAL::midpoint(hitsegmentCGAL.vertex(0), hitsegmentCGAL.vertex(1));
        hitPoint = Vector(CGAL::to_double(hitPointCGAL.x()),
                          CGAL::to_double(hitPointCGAL.y()),
                          CGAL::to_double(hitPointCGAL.z()));
      } else {
        throw GenerationErrorMessage(
            "This type of intersection should not happen");
      }

      LinkData& link = fluid->Links[iSolid];

      // This is set in any solid case
      float distanceInVoxels = (hitPoint - fluid->Position).GetMagnitude();
      // The distance is in voxels but must be output as a fraction of
      // the lattice vector. Scale it.
      link.Distance = distanceInVoxels / Neighbours::norms[iSolid];

      int ioletId = hitpoint_triangle_dist.first.second->id() - 2;
      // shifting back from unsigned.
      if (ioletId < 0) {
        // -1 => we hit a wall
        link.Type = geometry::CutType::WALL;
      } else {
        // We hit an inlet or outlet
        Iolet const& iolet = this->Iolets[ioletId];
        if (iolet.IsInlet) {
          link.Type = geometry::CutType::INLET;
        } else {
          link.Type = geometry::CutType::OUTLET;
        }
        // Set the Id
        link.IoletId = iolet.Id;
      }

      // If this link intersected the wall, store the normal of the cell we hit
      // and the distance to it.
      if (link.Type == geometry::CutType::WALL) {
        VectorCGAL CGALNorm = CGAL::cross_product(
            hitpoint_triangle_dist.first.second->halfedge()
                    ->next()
                    ->vertex()
                    ->point() -
                hitpoint_triangle_dist.first.second->halfedge()
                    ->vertex()
                    ->point(),
            hitpoint_triangle_dist.first.second->halfedge()
                    ->next()
                    ->next()
                    ->vertex()
                    ->point() -
                hitpoint_triangle_dist.first.second->halfedge()
                    ->next()
                    ->vertex()
                    ->point());
        CGALNorm = CGALNorm / CGAL::sqrt(CGALNorm.squared_length());
        link.WallNormalAtWallCut =
            Vector(CGALNorm.x(), CGALNorm.y(), CGALNorm.z());
        link.DistanceInVoxels = distanceInVoxels;
      }
    }
  }

  // If there's enough information available, an approximation of the wall
  // normal will be computed for this fluid site.
  this->ComputeAveragedNormal(site);
}

int PolyDataGenerator::Intersect(Site& site, Site& neigh) {
  int nHits;
  bool debugintersect = false;
  if (!neigh.IsFluidKnown) {
    // Neighbour unknown, must always intersect
    nHits = this->ComputeIntersectionsCGAL(site, neigh);
    if (nHits % 2 == 0) {
      // Even # hits, hence neigh has same type as site
      neigh.IsFluid = site.IsFluid;
    } else if (nHits % 2 == 1) {
      // Odd # hits, neigh is opposite type to site
      neigh.IsFluid = !site.IsFluid;
    } else {
      neigh.IsFluid = InsideOutside(neigh);
    }
    if (neigh.IsFluid)
      neigh.CreateLinksVector();
    neigh.IsFluidKnown = true;
  } else {
    // We know the fluidness of neigh, maybe don't need to intersect
    if (site.IsFluid != neigh.IsFluid) {
      nHits = this->ComputeIntersectionsCGAL(site, neigh);
      // Only in the case of difference must we intersect.
      if (nHits % 2 == 0) {
        bool Sinside = InsideOutside(site);
        bool Ninside = InsideOutside(neigh);
        throw InconsistentIntersectRayError(site, neigh, nHits, Sinside,
                                            Ninside);
      }
    } else {  // No need to intersect
      nHits = 0;
    }
  }

  if (debugintersect) {
    bool Sinside = InsideOutside(site);
    bool Ninside = InsideOutside(neigh);
    if (Sinside != site.IsFluid || Ninside != neigh.IsFluid) {
      throw InconsistentIntersectRayError(site, neigh, nHits, Sinside, Ninside);
    }
  }
  return nHits;
}

bool PolyDataGenerator::InsideOutside(Site& site) {
  PointCGAL point(site.Position[0], site.Position[1], site.Position[2]);
  bool inside;
  CGAL::Random_points_on_sphere_3<PointCGAL> random_point(1.);
  RayCGAL ray_query;
  int ori[3];
  bool nextray = true;
  int nHitsRay;
  PointCGAL v1;
  PointCGAL v2;
  PointCGAL v3;
  FacehandleCGAL f;
  std::vector<Object_and_primitive_id> rayhitcells;
  while (nextray) {
    ray_query = RayCGAL(point, *random_point);
    nHitsRay = this->AABBtree->number_of_intersected_primitives(ray_query);
    if (nHitsRay > 0) {
      rayhitcells.clear();
      this->AABBtree->all_intersections(ray_query,
                                        std::back_inserter(rayhitcells));
      for (std::vector<Object_and_primitive_id>::iterator i =
               rayhitcells.begin();
           i != rayhitcells.end(); ++i) {
        f = i->second;
        v1 = f->halfedge()->vertex()->point();
        v2 = f->halfedge()->next()->vertex()->point();
        v3 = f->halfedge()->next()->next()->vertex()->point();
        if (CGAL::orientation(v1, v2, v3, point) == 0) {
          nextray = false;
          inside = false;
          break;
        }
        ori[0] = CGAL::orientation(point, *random_point, v1, v2);
        ori[1] = CGAL::orientation(point, *random_point, v1, v3);
        ori[2] = CGAL::orientation(point, *random_point, v2, v3);
        if (ori[0] == 0 || ori[1] == 0 || ori[2] == 0) {
          nextray = true;
          ++random_point;
          break;
        }
        nextray = false;
        inside = nHitsRay % 2;
      }
    } else {
      inside = false;
      nextray = false;
    }
  }
  return inside;
}

int PolyDataGenerator::ComputeIntersections(Site& from, Site& to) {
  this->Locator->IntersectWithLine(&from.Position[0], &to.Position[0],
                                   this->hitPoints, this->hitCellIds);
  int hitpoints = this->hitPoints->GetNumberOfPoints();

  return hitpoints;
}

int PolyDataGenerator::ComputeIntersectionsCGAL(Site& from, Site& to) {
  PointCGAL p1(from.Position[0], from.Position[1], from.Position[2]);
  PointCGAL p2(to.Position[0], to.Position[1], to.Position[2]);
  PointCGAL p3;
  PointCGAL v1;
  PointCGAL v2;
  PointCGAL v3;
  FacehandleCGAL f;
  PointCGAL hitpoint;
  SegmentCGAL hitsegment;
  SegmentCGAL segment_query(p1, p2);
  int ori[5];
  int nHitsCGAL =
      this->AABBtree->number_of_intersected_primitives(segment_query);
  this->hitCellIdsCGAL.clear();
  this->IntersectionCGAL.clear();
  this->AABBtree->all_intersections(segment_query,
                                    std::back_inserter(this->hitCellIdsCGAL));
  Object_Primitive_and_distance OPD;

  if (nHitsCGAL) {
    for (std::vector<Object_and_primitive_id>::iterator i =
             this->hitCellIdsCGAL.begin();
         i != this->hitCellIdsCGAL.end(); ++i) {
      f = i->second;

      v1 = f->halfedge()->vertex()->point();
      v2 = f->halfedge()->next()->vertex()->point();
      v3 = f->halfedge()->next()->next()->vertex()->point();
      // 0,1,2 are not needed with only one point but the speedup looks minimal.
      ori[0] = CGAL::orientation(p1, p2, v1, v2);
      ori[1] = CGAL::orientation(p1, p2, v1, v3);
      ori[2] = CGAL::orientation(p1, p2, v2, v3);
      ori[3] = CGAL::orientation(p1, v1, v2, v3);
      ori[4] = CGAL::orientation(p2, v1, v2, v3);
      if (CGAL::assign(hitpoint, i->first)) {
        double distance =
            CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(hitpoint, p1)));
        OPD = std::make_pair(*i, distance);
        this->IntersectionCGAL.push_back(OPD);
      } else if (CGAL::assign(hitsegment, i->first)) {
        double distance1 = CGAL::to_double(
            CGAL::sqrt(CGAL::squared_distance(hitsegment.vertex(0), p1)));
        double distance2 = CGAL::to_double(
            CGAL::sqrt(CGAL::squared_distance(hitsegment.vertex(1), p1)));
        double distance = (distance1 + distance2) / 2;
        OPD = std::make_pair(*i, distance);
        this->IntersectionCGAL.push_back(OPD);
      } else {
        throw GenerationErrorMessage(
            "This type of intersection should not happen");
      }
      if (ori[0] == 0 || ori[1] == 0 || ori[2] == 0 || ori[3] == 0 ||
          ori[4] == 0) {
        // ori1,2,3 if the segment from voxel 1 to voxel 2 is in the same plane
        // as the edge. These 2 intersect and the result may be indetermined
        // Ori4 and Ori5. In this case either of the points are coplanar with
        // the triangle (primitive) In all cases the intersection is not
        // usefull.
        nHitsCGAL = -1;
      }
    }
  }

  if (nHitsCGAL != 1) {
    std::sort(this->IntersectionCGAL.begin(), this->IntersectionCGAL.end(),
              distancesort);
  }
  return nHitsCGAL;
}

// Function to be called on intersecting leaf nodes of the two OBB trees.
// Final void pointer is a pointer to an int, namely the count of the number
// of intersections found so far, which is incremented.
int IntersectingLeafCounter(vtkOBBNode* polyNode,
                            vtkOBBNode* cubeNode,
                            vtkMatrix4x4* transform,
                            void* ptr_to_intersection_count) {
  int& intersection_count = *static_cast<int*>(ptr_to_intersection_count);
  intersection_count++;
  return 0;
}

int PolyDataGenerator::BlockInsideOrOutsideSurface(const Block& block) {
  // Create an OBB tree for the block
  vtkSmartPointer<vtkOBBTree> blockSlightlyLargerOBBTree =
      block.CreateOBBTreeModel(1.0);

  // Count the number of domain OBB leaf nodes that intersect the single
  // node created for the block.
  int intersection_count = 0;
  Locator->IntersectWithOBBTree(blockSlightlyLargerOBBTree, NULL,
                                IntersectingLeafCounter,
                                static_cast<void*>(&intersection_count));
  // // Delete the underlying polydata
  // blockSlightlyLargerOBBTree->GetDataSet()->Delete();
  // // And the OBBTree itself
  // blockSlightlyLargerOBBTree->Delete();

  if (intersection_count == 0) {
    // either entirely inside or entirely outside
    double middlePosition[3];
    middlePosition[0] = block.Middle().Position[0];
    middlePosition[1] = block.Middle().Position[1];
    middlePosition[2] = block.Middle().Position[2];
    return Locator->InsideOrOutside(middlePosition);
  }
  return 0;
}

bool PolyDataGenerator::distancesort(const Object_Primitive_and_distance i,
                                     const Object_Primitive_and_distance j) {
  return (i.second < j.second);
}

}  // namespace hemelb::gmytool::gmy
