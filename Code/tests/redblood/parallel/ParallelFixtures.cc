// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/redblood/parallel/ParallelFixtures.h"

#include <random>

namespace hemelb::tests
{
    //! \brief gathers mid-domain and edge positions from all procs
    //! \details If there are insufficient number of edges, mid-domains are used instead.
    //! erase removes the components from the first process.
    std::vector<LatticeVector> GatherSpecialPositions(geometry::Domain const & domain,
						      size_t mid, size_t edges,
						      net::MpiCommunicator const &c)
    {
      std::random_device rd;
      std::mt19937 g(rd());

      int const nMids = domain.GetMidDomainSiteCount(0);
      int const nEdges = domain.GetDomainEdgeSiteCount(0);
      std::vector<LatticeVector> positions(c.Size() * (mid + edges));
      std::vector<int> shuf(nMids);
      std::iota(shuf.begin(), shuf.end(), 0);
      std::shuffle(shuf.begin(), shuf.end(), g);
      mid += std::max(0, static_cast<int>(edges) - nEdges);
      edges = std::min(edges, static_cast<size_t>(nEdges));
      for (size_t i(0); i < mid; ++i)
        {
          auto const site = domain.GetSite(shuf[i]);
          positions[c.Rank() * (mid + edges) + i] = site.GetGlobalSiteCoords();
        }
      shuf.resize(nEdges);
      std::iota(shuf.begin(), shuf.end(), 0);
      std::shuffle(shuf.begin(), shuf.end(), g);
      for (size_t i(0); i < edges; ++i)
        {
          auto const site = domain.GetSite(domain.GetMidDomainSiteCount() + shuf[i]);
          positions[c.Rank() * (mid + edges) + i + mid] = site.GetGlobalSiteCoords();
        }

      auto const sendType = net::MpiDataType<LatticePosition>();
      HEMELB_MPI_CALL(MPI_Allgather,
		      (MPI_IN_PLACE, mid + edges, sendType, positions.data(), mid + edges, sendType, c));

      // Erase contribution from first proc (acting as serial oracle)
      for (size_t i(0); i < positions.size() - mid - edges; ++i)
        {
          positions[i] = positions[i + mid + edges];
        }
      positions.resize(positions.size() - mid - edges);

      return positions;
    }

    DummyCell::DummyCell(LatticePosition const&position, LatticeForce f,
			 std::string const &templateName) :
      DummyCell(std::vector<LatticePosition> { position }, f, templateName)
    {
    }
    DummyCell::DummyCell(std::vector<LatticePosition> const &positions, LatticeForce f,
			 std::string const &templateName) :
      NodeCell(positions, templateName), force(f)
    {
    }

    LatticeEnergy DummyCell::operator()(std::vector<LatticeForceVector> &f) const
    {
      f.resize(GetNumberOfNodes());
      std::fill(f.begin(), f.end(), LatticeForceVector{force});
      return 0e0;
    }
    std::unique_ptr<redblood::CellBase> DummyCell::cloneImpl() const
    {
      return std::make_unique<DummyCell>(*this);
    }

    //! Creates list of cells for each set of positions from each process
    std::vector<hemelb::redblood::CellContainer::value_type> CreateCellsFromSpecialPositions(
											     geometry::Domain const & domain, size_t mid, size_t edges,
											     net::MpiCommunicator const &c, size_t nCells)
    {
      using hemelb::redblood::CellContainer;
      auto const positions = GatherSpecialPositions(domain, mid * nCells, edges * nCells, c);
      std::vector<CellContainer::value_type> cells;
      for (auto i_first = positions.begin(); i_first != positions.end(); i_first += mid + edges)
        {
          cells.push_back(std::make_shared<DummyCell>(std::vector<LatticePosition> { i_first,
		  i_first + mid
		  + edges },
	      1e0));
        }
      return cells;
    }

    net::MpiCommunicator CreateDumbGraphComm(net::MpiCommunicator const &comm)
    {
        // Dumb meaning that it communicates with every other rank
        auto neighbours = std::vector<int>(comm.Size() - 1);
        for (int i = 0, j = 0; i < comm.Size(); ++i) {
            if (i != comm.Rank())
                neighbours[j++] = i;
        }
        return comm.DistGraphAdjacent(neighbours);
    }
}
