// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include <deque>
#include <map>
#include <numeric>

#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>

#include "redblood/Mesh.h"
#include "Exception.h"
#include "constants.h"
#include "util/Iterator.h"

namespace hemelb::redblood
{
    // Check that the condition in Mesh.h is true.
    // TODO: fix as described in header
    static_assert(std::is_same<IdType, vtkIdType>::value,
		  "hemelb::redblood::IdType must be the same as vtkIdType");

    LatticePosition barycentre(MeshData::Vertices const &vertices)
    {
      typedef MeshData::Vertices::value_type Vertex;
      return std::accumulate(vertices.begin(), vertices.end(), Vertex(0, 0, 0))
          / Vertex::value_type(vertices.size());
    }
    LatticePosition barycentre(MeshData const &mesh)
    {
      return barycentre(mesh.vertices);
    }
    LatticeVolume volume(MeshData::Vertices const &vertices, MeshData::Facets const &facets)
    {
      auto result = std::transform_reduce(
              facets.begin(), facets.end(),
              LatticeVolume{0}, std::plus<LatticeVolume>{},
              [&](MeshData::Facet const& facet) {
                  auto& v0 = vertices[facet[0]];
                  auto& v1 = vertices[facet[1]];
                  auto& v2 = vertices[facet[2]];
                  return Dot(Cross(v0, v1), v2);
              }
      );
      // Minus sign comes from outward facing facet orientation
      return -result / 6.0;
    }
    LatticeVolume volume(MeshData const &mesh)
    {
      return volume(mesh.vertices, mesh.facets);
    }

    LatticeVolume area(MeshData::Vertices const &vertices, MeshData::Facets const &facets)
    {
        auto result = std::transform_reduce(
                facets.begin(), facets.end(),
                LatticeArea{0}, std::plus<LatticeArea>{},
                [&](MeshData::Facet const& facet) {
                    auto& v0 = vertices[facet[0]];
                    auto& v1 = vertices[facet[1]];
                    auto& v2 = vertices[facet[2]];
                    return Cross(v0 - v1, v2 - v1).GetMagnitude();
                }
        );
        return result * 0.5;
    }

    LatticeArea area(MeshData const &mesh)
    {
      return area(mesh.vertices, mesh.facets);
    }

    namespace
    {
      bool contains(MeshData::Facet const &a, MeshData::Facet::value_type v)
      {
        return a[0] == v or a[1] == v or a[2] == v;
      }
      bool edge_sharing(MeshData::Facet const &a, MeshData::Facet const &b)
      {
        return (contains(b, a[0]) ?
          1 :
          0) + (contains(b, a[1]) ?
          1 :
          0) + (contains(b, a[2]) ?
          1 :
          0) >= 2;
      }
      // Adds value as first non-negative number, if value not in array yet
      void insert(MeshData::Facet &container, MeshData::Facet::value_type value,
                  MeshData::Facet::value_type max)
      {
          for (auto& val: container) {
              if (val >= max)
              {
                  val = value;
                  return;
              }
              else if (val == value)
              {
                  return;
              }
          }
      }
    }

    MeshTopology::MeshTopology(MeshData const &mesh) :
            vertexToFacets(mesh.vertices.size()),
            facetNeighbors(mesh.facets.size())
    {
        // Loop over facets to create map from vertices to facets
        for (auto [i, facet]: util::enumerate(mesh.facets))
        {
            vertexToFacets.at(facet[0]).insert(i);
            vertexToFacets.at(facet[1]).insert(i);
            vertexToFacets.at(facet[2]).insert(i);
        }

        // Now creates map of neighboring facets
        IdType const N_FACETS = std::ssize(mesh.facets);
        std::array<IdType, 3> default_neigh = {N_FACETS, N_FACETS, N_FACETS};
        std::fill(facetNeighbors.begin(), facetNeighbors.end(), default_neigh);

        for (auto [i, facet]: util::enumerate(mesh.facets))
        {
            static_assert(std::is_same_v<decltype(facet), std::array<IdType,3> const&>);
            for (auto vertex_id: facet)
            {
                // check facets that this node is attached to
                auto const& neighboringFacets = vertexToFacets.at(vertex_id);
                for (auto neighboringFacet: neighboringFacets)
                {
                    if (i == neighboringFacet)
                        continue;

                    if (edge_sharing(facet, mesh.facets.at(neighboringFacet)))
                    {
                        insert(facetNeighbors.at(i), neighboringFacet, N_FACETS);
                    }
                }
            }
        }

#ifndef NDEBUG

      // Checks there are no uninitialized values
      for (auto const& facetNeighbor: facetNeighbors)
        for (unsigned int j = 0; j < 3; ++j)
        {
          assert(facetNeighbor[j] < N_FACETS);
        }

#endif
    }

    void Mesh::operator*=(Dimensionless const &scale)
    {
      auto const barycentre = GetBarycentre();

      for (auto &vertex : mesh->vertices)
      {
        vertex = (vertex - barycentre) * scale + barycentre;
      }
    }

    void Mesh::operator*=(util::Matrix3D const &rotation)
    {
      auto const barycentre = GetBarycentre();

      for (auto &vertex : mesh->vertices)
      {
        rotation.timesVector(vertex - barycentre, vertex);
        vertex += barycentre;
      }
    }

    void Mesh::operator+=(LatticePosition const &offset)
    {
      for (auto &vertex : mesh->vertices)
      {
        vertex += offset;
      }
    }

    void Mesh::operator+=(std::vector<LatticePosition> const &displacements)
    {
      assert(displacements.size() == mesh->vertices.size());
      auto i_disp = displacements.begin();

      for (auto &vertex : mesh->vertices)
      {
        vertex += * (i_disp++);
      }
    }

    namespace
    {
      std::shared_ptr<MeshData> initial_tetrahedron()
      {
        std::shared_ptr<MeshData> data(new MeshData);

        // facets at something degrees from one another
        data->vertices.push_back(LatticePosition(0, 0, 0));
        data->vertices.push_back(LatticePosition(1, 0, 1));
        data->vertices.push_back(LatticePosition(1, 1, 0));
        data->vertices.push_back(LatticePosition(0, 1, 1));

        redblood::MeshData::Facet indices;
        indices[0] = 0;
        indices[1] = 2;
        indices[2] = 3;
        data->facets.push_back(indices);
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        data->facets.push_back(indices);
        indices[0] = 0;
        indices[1] = 3;
        indices[2] = 1;
        data->facets.push_back(indices);
        indices[0] = 1;
        indices[1] = 3;
        indices[2] = 2;
        data->facets.push_back(indices);
        return data;
      }

      size_t vertex(std::shared_ptr<MeshData> &data,
                    std::map<std::pair<size_t, size_t>, size_t> &vertices, size_t const &i0,
                    size_t const &i1)
      {
        std::pair<size_t, size_t> const indices(i0, i1);
        std::map<std::pair<size_t, size_t>, size_t>::const_iterator i_found =
            vertices.find(indices);

        if (i_found == vertices.end())
        {
          i_found = vertices.find(std::pair<size_t, size_t>(i1, i0));
        }

        if (i_found == vertices.end())
        {
          data->vertices.push_back( (data->vertices[i0] + data->vertices[i1]) * 0.5);
          vertices[indices] = data->vertices.size() - 1;
          return data->vertices.size() - 1;
        }

        return i_found->second;
      }

      void refine(std::shared_ptr<MeshData> &data)
      {
        MeshData::Facets const facets(data->facets);
        data->facets.clear();
        data->facets.resize(facets.size() * 4);
        data->vertices.reserve(data->facets.size() * 3 + data->vertices.size());

        // Container with midpoint indices, so midpoints are only added once
        typedef std::pair<size_t, size_t> Pair;
        std::map<Pair, size_t> new_vertices;

        MeshData::Facets::const_iterator i_orig_facet(facets.begin());
        MeshData::Facets::const_iterator const i_orig_facet_end(facets.end());
        MeshData::Facets::iterator i_facet = data->facets.begin();

        for (; i_orig_facet != i_orig_facet_end; ++i_orig_facet)
        {
          MeshData::Facet::value_type const i0 = (*i_orig_facet)[0];
          MeshData::Facet::value_type const i1 = (*i_orig_facet)[1];
          MeshData::Facet::value_type const i2 = (*i_orig_facet)[2];

          // Adds new vertices halfway through edges
          MeshData::Facet::value_type mid0(vertex(data, new_vertices, i0, i1)),
              mid1(vertex(data, new_vertices, i1, i2)), mid2(vertex(data, new_vertices, i2, i0));

          // Adds all four new faces
          (*i_facet)[0] = i0;
          (*i_facet)[1] = mid0;
          (*i_facet)[2] = mid2;
          ++i_facet;

          (*i_facet)[0] = mid0;
          (*i_facet)[1] = i1;
          (*i_facet)[2] = mid1;
          ++i_facet;

          (*i_facet)[0] = mid1;
          (*i_facet)[1] = i2;
          (*i_facet)[2] = mid2;
          ++i_facet;

          (*i_facet)[0] = mid0;
          (*i_facet)[1] = mid1;
          (*i_facet)[2] = mid2;
          ++i_facet;
        }
      }
    }

    Mesh refine(Mesh data, unsigned int depth)
    {
      if (depth == 0)
      {
        return data.clone();
      }

      std::shared_ptr<MeshData> newData(new MeshData(*data.GetData()));

      for (unsigned int i(0); i < depth; ++i)
      {
        refine(newData);
      }

      return Mesh(newData);
    }

    Mesh tetrahedron(unsigned int depth)
    {
      std::shared_ptr<MeshData> result(initial_tetrahedron());

      for (unsigned int i(0); i < depth; ++i)
      {
        refine(result);
      }

      return Mesh(result);
    }

    Mesh pancakeSamosa(unsigned int depth)
    {
      std::shared_ptr<redblood::MeshData> mesh(new redblood::MeshData);

      // facets at something degrees from one another
      mesh->vertices.push_back(LatticePosition(0, 0, 0));
      mesh->vertices.push_back(LatticePosition(1, 0, 1));
      mesh->vertices.push_back(LatticePosition(1, 1, 0));

      redblood::MeshData::Facet indices;
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      mesh->facets.push_back(indices);
      indices[0] = 2;
      indices[1] = 1;
      indices[2] = 0;
      mesh->facets.push_back(indices);

      // Create topology by hand cos we generally don't allow for this kind of
      // ambiguous self-referencing shape.
      std::shared_ptr<redblood::MeshTopology> topo(new redblood::MeshTopology);
      MeshTopology::VertexToFacets::value_type v2f;
      v2f.insert(0);
      v2f.insert(1);
      topo->vertexToFacets.resize(3, v2f);

      MeshTopology::FacetNeighbors::value_type neighbors[2] = { { { 0, 0, 0 } }, { { 1, 1, 1 } } };
      topo->facetNeighbors.push_back(neighbors[1]);
      topo->facetNeighbors.push_back(neighbors[0]);

      return refine(Mesh(mesh, topo), depth);
    }

    Mesh icoSphere(unsigned int depth)
    {
      std::shared_ptr<redblood::MeshData> mesh(new redblood::MeshData);
      // First creates the simplest icosphere
      // Dumbly copied from
      // http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
      LatticeDistance const t = 0.5 * (1. + std::sqrt(5.0));
      mesh->vertices.clear();
      mesh->vertices.emplace_back(-1, t, 0);
      mesh->vertices.emplace_back(1, t, 0);
      mesh->vertices.emplace_back(-1, -t, 0);
      mesh->vertices.emplace_back(1, -t, 0);
      mesh->vertices.emplace_back(0, -1, t);
      mesh->vertices.emplace_back(0, 1, t);
      mesh->vertices.emplace_back(0, -1, -t);
      mesh->vertices.emplace_back(0, 1, -t);
      mesh->vertices.emplace_back(t, 0, -1);
      mesh->vertices.emplace_back(t, 0, 1);
      mesh->vertices.emplace_back(-t, 0, -1);
      mesh->vertices.emplace_back(-t, 0, 1);
      auto intel_compiler = [](size_t i, size_t j, size_t k) -> MeshData::Facet
      {
        MeshData::Facet result;
        result[0] = i; result[1] = j; result[2] = k;
        return result;
      };
      mesh->facets.clear();
      mesh->facets.push_back(intel_compiler(0, 11, 5));
      mesh->facets.push_back(intel_compiler(0, 5, 1));
      mesh->facets.push_back(intel_compiler(0, 1, 7));
      mesh->facets.push_back(intel_compiler(0, 7, 10));
      mesh->facets.push_back(intel_compiler(0, 10, 11));
      mesh->facets.push_back(intel_compiler(1, 5, 9));
      mesh->facets.push_back(intel_compiler(5, 11, 4));
      mesh->facets.push_back(intel_compiler(11, 10, 2));
      mesh->facets.push_back(intel_compiler(10, 7, 6));
      mesh->facets.push_back(intel_compiler(7, 1, 8));
      mesh->facets.push_back(intel_compiler(3, 9, 4));
      mesh->facets.push_back(intel_compiler(3, 4, 2));
      mesh->facets.push_back(intel_compiler(3, 2, 6));
      mesh->facets.push_back(intel_compiler(3, 6, 8));
      mesh->facets.push_back(intel_compiler(3, 8, 9));
      mesh->facets.push_back(intel_compiler(4, 9, 5));
      mesh->facets.push_back(intel_compiler(2, 4, 11));
      mesh->facets.push_back(intel_compiler(6, 2, 10));
      mesh->facets.push_back(intel_compiler(8, 6, 7));
      mesh->facets.push_back(intel_compiler(9, 8, 1));
      orientFacets(*mesh);
      // Then refines it
      for (unsigned int i(0); i < depth; ++i)
      {
        refine(mesh);
      }
      // Finally, makes sure vertices are on unit-sphere
      for (auto &vertex : mesh->vertices)
      {
        vertex.Normalise();
      }
      return mesh;
    }

    unsigned orientFacets(Mesh &mesh, bool outward)
    {
      return orientFacets(*mesh.GetData(), outward);
    }

    using UnitVector = util::Vector3D<Dimensionless>;

    // What index (i.e. 0, 1, 2) does the point `id` have in the facet `f`?
    // If not present will return 3.
    // This is known as Facet = std::array<IdType, 3>;
    unsigned find_index(MeshData::Facet const& f, IdType const& id) {
      return (unsigned)std::distance(f.begin(), std::find(f.begin(), f.end(), id));
    }

    struct MeshConnectivity {
      IdType nPts;
      IdType nTris;
      
      // We represent point-to-cell data with two arrays of `IdType`. First:
      //   std::vector<IdType> offsets;
      // which will hold at index `i` the offset into the second array
      // where the data for point `i` can be found.
      //
      // Second: this array holds the concatenated data for all points.
      // One point's data consists of:
      //  - first, the number of triangles that use it
      //  - second, the ids of all triangles that use it
      // We know that there are 3*nTri uses of points + 1 value per
      // point to hold the count.
      //   std::vector<IdType> tris_using_pt;
      //
      // But, we can stick these together offsets first, followed by tris_using_pt
      std::vector<IdType> point2cell;

      // We have triangles on a closed, manifold surface, so each
      // triangle connects to exactly three others.
      std::vector<IdType> neighbours;

      IdType& neigh(IdType triId, IdType i) {
	return neighbours[3*triId + i];
      }

      IdType const& neigh(IdType triId, IdType i) const {
	return neighbours[3*triId + i];
      }
      
      explicit MeshConnectivity(MeshData const& mesh) :
	nPts(mesh.vertices.size()),
	nTris(mesh.facets.size()),
	point2cell(2*nPts + 3*nTris),
	neighbours(3*nTris, -1)
	// point_starts(nPts),
	// point_to_cell(nTris * 3 + nPts)
      {
	// First get the use counts of each point
	// 
	// Initialise in the IIFE so can have const
	auto const point_use_counts = [&]() {
	  auto ans = std::vector<IdType>(nPts, 0);
	  for (auto const& tri: mesh.facets) {
	    for (int dim = 0; dim < 3; ++dim) {
	      ++ans[tri[dim]];
	    }
	  }
	  return ans;
	}();

	// Compute the scan/partial sum to get the offsets (if we
	// didn't have counts). We really want the exclusive_scan, but
	// don't have that until C++17
	auto const cumulative_use_count = [&]() {
	  auto ans = std::vector<IdType>(nPts);
	  ans[0] = 0;
	  std::partial_sum(point_use_counts.begin(), --point_use_counts.end(),
			   ++ans.begin());
	  return ans;
	}();
	// Know sum(use_counts) == nTris * 3 and last elem of
	// cumulative use counts is sum of all except the last
	// element.
	if (cumulative_use_count[nPts - 1] + point_use_counts[nPts - 1]
	    != nTris * 3) {
	  throw Exception() << "Use counts incorrect. I can't use partial_sum right.";
	}

	// Initialise point_starts, accounting for per-point count,
	// and zero the count in output (as we shall use that to
	// figure out where to write the cell IDs as we go).
	for (IdType ptId = 0; ptId < nPts; ++ptId) {
	  point2cell[ptId] = cumulative_use_count[ptId] + ptId + nPts;
	  point2cell[point2cell[ptId]] = 0;
	}

	// Now fill point2cell and neighbours
	for (IdType triId = 0; triId < nTris; ++triId) {
	  for (int dim = 0; dim < 3; ++dim) {
	    // For each point in each triangle, get its ID
	    auto const& ptId = mesh.facets[triId][dim];

	    // Where's its data?
	    auto const& start = point2cell[ptId];
	    // Note mutable ref for inc below
	    auto& nSeen = point2cell[start];

	    // Edge i goes from point i -> i+1
	    auto const& endPtId = mesh.facets[triId][(dim + 1) % 3];

	    // First loop over tris we already know share ptId with
	    // us, looking for those that also share endPtId. Those
	    // are edge-neighbours. We do this first to avoid having
	    // to exclude ourself.
	    for (IdType const& neighId: iter_tris_using_pt(ptId)) {
	      const auto& neighFacet = mesh.facets[neighId];
	      auto neighEndDim = find_index(neighFacet, endPtId);
	      if (neighEndDim != neighFacet.size()) {
		// Found, so this is an edge neighbour, record as such
		if (neigh(triId, dim) != -1)
		  throw Exception() << "Setting a neighbour map value that has already been set";

		neigh(triId, dim) = neighId;
		// Find the other end within the neighbouring tri to
		// figure out which way round it is.
		auto neighDim = find_index(neighFacet, ptId);
		if (neighDim == neighFacet.size())
		  throw Exception() << "Can't find a point that I really should";

		auto neighEdgeId = ((neighDim + 1) % 3 == neighEndDim) ?
		  // This is edge neighDim
		  neighDim :
		  // This is edge neighDimEnd
		  neighEndDim;

		if (neigh(neighId, neighEdgeId) != -1)
		  throw Exception() << "Setting a neighbour map value that has already been set";
		neigh(neighId, neighEdgeId) = triId;
	      }
	    }

	    point2cell[start + 1 + nSeen] = triId;
	    ++nSeen;

	  }
	}

#ifndef NDEBUG
	// Check counts are correct
	for (IdType ptId = 0; ptId < nPts; ++ptId) {
	  if (point2cell[point2cell[ptId]] != point_use_counts[ptId]) {
	    throw Exception() << "Point use count not correct";
	  }
	}

	// Check no neighbours uninitialised
	for (IdType triId = 0; triId < nTris; ++triId) {
	  for (int i = 0; i < 3; ++i) {
	    if (neighbours[3*triId + i] == -1) {
	      throw Exception()
		<< "Uninitialised neighbour map for triangle " << triId << " edge " << i;
	    }
	  }
	}
#endif

      }

      struct pt_cell_range {
	IdType const * _begin;
	IdType const * _end;

	pt_cell_range(MeshConnectivity const* mc, IdType ptId) {
	  auto data_idx = mc->point2cell[ptId];
	  auto n = mc->point2cell[data_idx];

	  _begin = &mc->point2cell[data_idx + 1];
	  _end = &mc->point2cell[data_idx + 1 + n];
	}
	IdType const * begin() const {
	  return _begin;
	}
	IdType const* end() const {
	  return _end;
	}
      };

      pt_cell_range iter_tris_using_pt(IdType ptId) const {
	return pt_cell_range{this, ptId};
      }

    };

    unsigned orientFacets(MeshData &mesh, bool outward)
    {
      auto const nTris = mesh.facets.size();

      // Mutable as we'll be flipping facets below
      MeshConnectivity conn{mesh};

      std::vector<UnitVector> normals(nTris);
      std::transform(mesh.facets.cbegin(), mesh.facets.cend(),
		     normals.begin(),
		     [&mesh](MeshData::Facet const& facet) {
		       auto& v0 = mesh.vertices[facet[0]];
		       auto& v1 = mesh.vertices[facet[1]];
		       auto& v2 = mesh.vertices[facet[2]];
               return Cross(v0 - v1, v2 - v1).GetNormalised();
		     });

      // Do as vtkPolyDataNormals, does, but simplified as only have
      // one connected component.
      // 
      // The left-most polygon should have its outward pointing normal
      // facing left. If it doesn't, reverse the vertex order. Then
      // use it as the seed for other connected polys.

      // To find left-most polygon, first find left-most point, and
      // examine neighboring polys and see which one has a normal
      // that's "most aligned" with the X-axis.
      auto leftmost_pt_iter =
	std::min_element(mesh.vertices.cbegin(), mesh.vertices.cend(),
			 [](LatticePosition const& a, LatticePosition const& b) {
			   return a.x() < b.x();
			 });
      const IdType leftmost_pt_id = std::distance(mesh.vertices.cbegin(), leftmost_pt_iter);
      auto leftmost_cells = conn.iter_tris_using_pt(leftmost_pt_id);
      if (leftmost_cells.begin() == leftmost_cells.end()) {
	throw Exception() << "No leftmost triangles found when orienting facets";
      }
	
      auto leftmost_tri_id =
	*std::max_element(leftmost_cells.begin(), leftmost_cells.end(),
			  [&normals](std::size_t a, std::size_t b) {
			    return std::fabs(normals[a].x()) < std::fabs(normals[b].x());
			  });

      unsigned nFlips = 0;
      auto flip = [&](IdType triId) {
	// Half-edge i is from pt i -> (i+1)%3
	// 
	//  0        2
	//  |\ 2     |\ 2
	//  | \      | \
	// 0|  2 => 1|  0
	//  | /      | /
	//  |/ 1     |/ 0
	//  1        1
	std::swap(mesh.facets[triId][0], mesh.facets[triId][2]);
	std::swap(conn.neigh(triId,0), conn.neigh(triId, 1));
	normals[triId] = -normals[triId];
	nFlips++;
      };

      // We need to track whether cells have been checked
      auto visited = std::vector<bool>(nTris, false);

      // Flip if pointing to the right.
      if (normals[leftmost_tri_id].x() > 0) {
	flip(leftmost_tri_id);
      }
      visited[leftmost_tri_id] = true;

      // Now we have at least one properly oriented triangle, spread this.       
      auto edge_tris = std::deque<IdType>{};
      edge_tris.push_back(leftmost_tri_id);

      while (edge_tris.size()) {
	IdType const curTriId = edge_tris.front();
	edge_tris.pop_front();

	auto& curTri = mesh.facets[curTriId];

	// For each edge on the triangle
	for (int i1 = 0; i1 < 3; i1++) {
	  int i2 = (i1 + 1) % 3;
	  IdType p1 = curTri[i1];
	  IdType p2 = curTri[i2];

	  // Find the triangle that shares this edge.
	  IdType const& neighId = conn.neigh(curTriId, i1);

	  if (visited[neighId])
	    continue;
	  auto const& neigh = mesh.facets[neighId];

	  // What order are the points in?
	  int j1 = find_index(neigh, p1);
	  int j2 = find_index(neigh, p2);
	  // If the edges run the same way, then neigh needs to flip,
	  // otherwise fine. By that we mean:
	  if ((j1 + 1) % 3 == j2) {
	    flip(neighId);
	  }
	  visited[neighId] = true;
	  edge_tris.push_back(neighId);
	}
      }
#ifndef NDEBUG
      if (!std::all_of(visited.begin(), visited.end(), [](bool x) { return x; }))
	throw Exception() << "Some triangles not visited";
#endif

      return nFlips;
    }

    unsigned orientFacets(MeshData &mesh, vtkPolyData &polydata, bool outward)
    {
      vtkSmartPointer<vtkPolyDataNormals> normalCalculator = vtkSmartPointer<vtkPolyDataNormals>::New();
      normalCalculator->SetInputData(&polydata);
      normalCalculator->ComputePointNormalsOff();
      normalCalculator->ComputeCellNormalsOn();
      normalCalculator->SetAutoOrientNormals(1);
      normalCalculator->Update();
      auto normals = normalCalculator->GetOutput()->GetCellData()->GetNormals();

      if (normals->GetNumberOfComponents() != 3)
	throw Exception() << "Normals does not have 3 components";

      if (normals->GetNumberOfTuples() != IdType(mesh.facets.size()))
	throw Exception() << "Normals does not have " << mesh.facets.size() << " entries, it has " << normals->GetNumberOfTuples();

      // Loop over each facet, checks orientation and modify as appropriate
      unsigned normalId = 0;
      unsigned numSwapped = 0;
      for (auto &facet : mesh.facets)
      {
        auto const &v0 = mesh.vertices[facet[0]];
        auto const &v1 = mesh.vertices[facet[1]];
        auto const &v2 = mesh.vertices[facet[2]];

        double vtkNormal[3];
        normals->GetTuple(normalId, vtkNormal);
        LatticePosition direction(vtkNormal[0], vtkNormal[1], vtkNormal[2]);

        if ( ( Dot(Cross(v0 - v1, v2 - v1), direction) > 0e0) xor outward)
        {
          std::swap(facet[0], facet[2]);
          ++numSwapped;
        }

        ++normalId;
      }

      return numSwapped;
    }

} // hemelb::redblood
