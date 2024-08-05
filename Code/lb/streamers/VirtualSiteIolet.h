// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_VIRTUALSITEIOLET_H
#define HEMELB_LB_STREAMERS_VIRTUALSITEIOLET_H

#include <map>
#include <fstream>

#include <boost/container/flat_map.hpp>

#include "geometry/neighbouring/RequiredSiteInformation.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "lb/lattices/LatticeInfo.h"
#include "lb/streamers/VirtualSite.h"
#include "log/Logger.h"

namespace hemelb::lb
{

    template<collision_type C>
    class VirtualSiteIolet
    {
    public:
        using CollisionType = C;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename CollisionType::LatticeType;

    private:
        using VSiteType = VirtualSite<LatticeType>;

          CollisionType collider;
          BulkLink<CollisionType> bulkLinkDelegate;
          BounceBackLink<CollisionType> wallLinkDelegate;
          BoundaryValues* bValues;
          const geometry::neighbouring::NeighbouringDomain& neighbouringLatticeData;

          // These will store a map from localIdx => (iolet, vsite, direction) triples
          struct IoletVSiteDirection
          {
              IoletVSiteDirection(InOutLet*iolet_, VirtualSite<LatticeType>* vsite_, Direction i_) :
                  iolet(iolet_), vsite(vsite_), direction(i_)
              {
              }
              InOutLet* iolet;
              VirtualSite<LatticeType>* vsite;
              Direction direction;
          };
          using VSiteByLocalIdxMultiMap = boost::container::flat_multimap<site_t, IoletVSiteDirection>;
          VSiteByLocalIdxMultiMap vsByLocalIdx;

        public:
          VirtualSiteIolet(InitParams& initParams) :
              collider(initParams), bulkLinkDelegate(collider, initParams),
                  wallLinkDelegate(collider, initParams), bValues(initParams.boundaryObject),
                  neighbouringLatticeData(initParams.latDat->GetNeighbouringData())
          {
            // Loop over the local in/outlets, creating the extra data objects.
            unsigned nIolets = bValues->GetLocalIoletCount();
            for (unsigned iIolet = 0; iIolet < nIolets; ++iIolet)
            {
              InOutLet& iolet = *bValues->GetLocalIolet(iIolet);
              if (iolet.GetExtraData() == nullptr)
                iolet.SetExtraData(new VSExtra<LatticeType>(iolet));
            }

            auto& lattice = LatticeType::GetLatticeInfo();
            // Want to loop over each site this streamer is responsible for,
            // as specified in the siteRanges.
            for (auto [start, end]: initParams.siteRanges)
            {
              for (site_t siteIdx = start; siteIdx < end; ++siteIdx)
              {
                geometry::Site<const geometry::Domain> site =
                    initParams.latDat->GetSite(siteIdx);

                if (site.GetSiteType() != bValues->GetIoletType())
                {
                  std::cout << "  Site is not of correct iolet type " << siteIdx << std::endl;
                  continue;
                }
                const LatticeVector siteLocation = site.GetGlobalSiteCoords();

                // Get the iolet
                InOutLet& iolet = *bValues->GetLocalIolet(site.GetIoletId());

                // Get the extra data for this iolet
                VSExtra<LatticeType>* extra = GetExtra(&iolet);

                // For each site index in the neighbourhood
                for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i)
                {
                  // If this direction doesn't cross an iolet, skip it.
                  if (!site.HasIolet(i))
                    continue;

                  const LatticeVector neighbourLocation = siteLocation + lattice.GetVector(i);
                  site_t neighbourGlobalIdx =
                      initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(neighbourLocation);
                  typename VSiteType::Map::iterator vNeigh = extra->vSites.find(neighbourGlobalIdx);

                  // find() returns end() if key not present
                  if (vNeigh == extra->vSites.end())
                  {
                    // Create a vSite
                    std::pair<typename VSiteType::Map::iterator, bool> inserted =
                        extra->vSites.insert(typename VSiteType::Map::value_type(neighbourGlobalIdx,
                                                                                 VSiteType(initParams,
                                                                                           *extra,
                                                                                           neighbourLocation)));
                    // inserted.first is an iterator, pointing to a pair<neighGlobalIdx, newly constructed vSite>
                    vNeigh = inserted.first;
                  }

                  // Add the (possibly newly created) virtual site to the map by local index.
                  vsByLocalIdx.insert(typename VSiteByLocalIdxMultiMap::value_type(siteIdx,
                                                                                   IoletVSiteDirection(&iolet,
                                                                                                       & (vNeigh->second),
                                                                                                       lattice.GetInverseIndex(i))));
                }
              }
            }

          }

          /*
           * Stream and collide as normal for bulk and wall links. The iolet
           * links will be done in the post-step as we must ensure that all
           * the data is available to construct virtual sites.
           */
          void StreamAndCollide(const site_t firstIndex, const site_t lastIndex,
                                         const LbmParameters* lbmParams,
                                         geometry::FieldData& latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            auto&& dom = latDat.GetDomain();
            for (site_t siteIndex = firstIndex; siteIndex < lastIndex; siteIndex++)
            {
              auto&& site = latDat.GetSite(siteIndex);
              VarsType hydroVars(site);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
              {
                if (site.HasWall(ii))
                {
                  wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
                else if (site.HasIolet(ii))
                {
                  // Do nothing for now.
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }
              }
              /*
               * Store the density and velocity for later use.
               */
              // Get the iolet
              InOutLet* iolet = bValues->GetLocalIolet(site.GetIoletId());

              // Get the extra data for this iolet
              VSExtra<LatticeType>* extra = GetExtra(iolet);
              site_t globalIdx = dom.GetGlobalNoncontiguousSiteIdFromGlobalCoords(site.GetGlobalSiteCoords());
              RSHV& cachedHV = extra->hydroVarsCache[globalIdx];
              cachedHV.t = bValues->GetTimeStep();
              cachedHV.rho = hydroVars.density;
              cachedHV.u = hydroVars.velocity;

                UpdateCachePostCollision(site,
                                         hydroVars,
                                         lbmParams,
                                         propertyCache);
            }
          }

          void PostStep(const site_t firstIndex, const site_t lastIndex,
                                 const LbmParameters* lbmParams, geometry::FieldData* latDat,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            const LatticeTimeStep t = bValues->GetTimeStep();
            const typename VSiteByLocalIdxMultiMap::iterator beginVSites =
                vsByLocalIdx.lower_bound(firstIndex), endVSites =
                vsByLocalIdx.lower_bound(lastIndex);

            for (typename VSiteByLocalIdxMultiMap::iterator vSiteIt = beginVSites;
                vSiteIt != endVSites; ++vSiteIt)
            {
              site_t siteIdx = vSiteIt->first;
              // vSiteIt->second == (Iolet*, VirtualSite*, Direction)
              InOutLet* iolet = vSiteIt->second.iolet;
              // Get the extra data for this iolet
              VSExtra<LatticeType>* extra = GetExtra(iolet);
              VSiteType* vSite = vSiteIt->second.vsite;

              // Compute the distributions for the vSite if needed
              CalculateVirtualSiteDistributions(*latDat, *iolet, extra->hydroVarsCache, *vSite, t);
              // Stream this direction
              Direction i = vSiteIt->second.direction;
              * (latDat->GetFNew(siteIdx * LatticeType::NUMVECTORS + i)) = vSite->hv.fPostColl[i];
              //* (m_fieldData->GetFNew(GetBBIndex(site.GetIndex(), direction))) = hydroVars.GetFPostCollision()[direction];
              //return (siteIndex * LatticeType::NUMVECTORS) + LatticeType::INVERSEDIRECTIONS[direction];
            }
          }

          static void DumpTables(const VirtualSiteIolet* ioletStreamer,
                                 const VirtualSiteIolet* ioletWallStreamer,
                                 const geometry::Domain* latDat)
          {
            InOutLet* iolet = ioletStreamer->bValues->GetLocalIolet(0);
            VSExtra<LatticeType>* extra = GetExtra(iolet);

            std::ofstream hvCache("hvCache");
            hvCache << "# local global x y z" << std::endl;
            for (RSHV::Map::const_iterator hvIt = extra->hydroVarsCache.begin();
                hvIt != extra->hydroVarsCache.end(); ++hvIt)
            {
              site_t global = hvIt->first;
              LatticeVector pos;
              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(global, pos);
              site_t local = latDat->GetContiguousSiteId(pos);
              hvCache << local << " " << global << " " << pos << std::endl;
            }
            hvCache.close();

            std::ofstream vSites("vSites");
            vSites << "# global x y z vSitePtr" << std::endl;
            for (typename VirtualSite<LatticeType>::Map::const_iterator vsIt =
                extra->vSites.begin(); vsIt != extra->vSites.end(); ++vsIt)
            {
              site_t global = vsIt->first;
              LatticeVector pos;
              const VSiteType& vs = vsIt->second;

              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(global, pos);
              vSites << global << " " << pos << " " << &vs;
              for (unsigned i = 0; i < vs.neighbourGlobalIds.size(); ++i)
              {
                vSites << " " << vs.neighbourGlobalIds[i];
              }
              vSites << std::endl;
            }
            vSites.close();

            std::ofstream outletMap("outletMap");
            outletMap << "# local global x y z vSitePtr direction" << std::endl;
            for (typename VSiteByLocalIdxMultiMap::const_iterator entry =
                ioletStreamer->vsByLocalIdx.begin(); entry != ioletStreamer->vsByLocalIdx.end();
                ++entry)
            {
              site_t local = entry->first;
              geometry::Site<const geometry::Domain> site = latDat->GetSite(local);
              LatticeVector pos = site.GetGlobalSiteCoords();
              site_t global = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(pos);
              outletMap << local << " " << global << " " << pos
                  << " " << entry->second.vsite << " " << entry->second.direction << std::endl;
            }
            outletMap.close();

            std::ofstream outletWallMap("outletWallMap");
            outletWallMap << "# local global x y z vSitePtr direction" << std::endl;
            for (typename VSiteByLocalIdxMultiMap::const_iterator entry =
                ioletWallStreamer->vsByLocalIdx.begin();
                entry != ioletWallStreamer->vsByLocalIdx.end(); ++entry)
            {
              site_t local = entry->first;
              geometry::Site<const geometry::Domain> site = latDat->GetSite(local);
              LatticeVector pos = site.GetGlobalSiteCoords();
              site_t global = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(pos);
              outletWallMap << local << " " << global << " " << pos << " " << entry->second.vsite << " " << entry->second.direction
                  << std::endl;
            }
          }

        private:
          static VSExtra<LatticeType>* GetExtra(InOutLet* iolet)
          {
            // Get the extra data for this iolet
            auto* ans = dynamic_cast<VSExtra<LatticeType>*>(iolet->GetExtraData());
            if (ans == nullptr)
            {
              // panic
              log::Logger::Log<log::Critical, log::OnePerCore>("Extra data not available for in/outlet. Aborting.\n");
              std::exit(1);
            }
            return ans;
          }

          void CalculateVirtualSiteDistributions(const geometry::FieldData& latDat,
                                                 const InOutLet& iolet, RSHV::Map& hydroVarsCache,
                                                 VSiteType& vSite, const LatticeTimeStep t)
          {
            if (vSite.hv.t != t)
            {
              vSite.hv.rho = CalculateVirtualSiteDensity(latDat, iolet, hydroVarsCache, vSite, t);
              vSite.hv.u = CalculateVirtualSiteVelocity(latDat, iolet, hydroVarsCache, vSite, t);
              vSite.hv.t = t;

              // Should really compute stress, relax it with collision and
              // then stream that, but will just use the equilibrium for now.
              // distribn_t fEq[LatticeType::NUMVECTORS];
              LatticeType::CalculateFeq(vSite.hv.rho,
                                        vSite.hv.u,
                                        vSite.hv.fPostColl);
            }
          }

          /**
           * Choose the density of the virtual site such that the
           * interpolated error (in a sum of squares sense) at the iolet
           * plane is minimised.
           *
           * rho_virtual = SUM_i ( q[i] (rho_iolet - (1 - q[i]) * rho[i]) ) / SUM _i (q[i]^2)
           *
           * @param latDat
           * @param iolet
           * @param hydroVarsCache
           * @param vSite
           * @param t
           * @return
           */
          LatticeDensity CalculateVirtualSiteDensity(const geometry::FieldData& latDat,
                                                     const InOutLet& iolet,
                                                     RSHV::Map& hydroVarsCache,
                                                     const VSiteType& vSite,
                                                     const LatticeTimeStep t)
          {
            LatticeDensity rho = 0.;
            LatticeDensity rho_iolet = iolet.GetDensity(t);
            for (unsigned i = 0; i < vSite.neighbourGlobalIds.size(); ++i)
            {

              LatticeDensity rho_site_i = GetHV(latDat,
                                                hydroVarsCache,
                                                vSite.neighbourGlobalIds[i],
                                                t).rho;
              rho += vSite.q[i] * (rho_iolet - (1.0 - vSite.q[i]) * rho_site_i);
            }
            rho /= vSite.sumQiSq;
            return rho;
          }
          /**
           * Choose the velocity given by a least squares fit of the neighbouring
           * points' normal velocities. Normal only, since the velocity beyond
           * the iolet is assumed to be laminar Poiseuille flow.
           *
           * Anstaz is u(x,y) = Ax + By + C
           *
           * @param latDat
           * @param iolet
           * @param hydroVarsCache
           * @param vSite
           * @param t
           * @return
           */
          LatticeVelocity CalculateVirtualSiteVelocity(const geometry::FieldData& latDat,
                                                       const InOutLet& iolet,
                                                       RSHV::Map& hydroVarsCache,
                                                       const VSiteType& vSite,
                                                       const LatticeTimeStep t)
          {
            /*
             *  Anstaz is u(x,y) = Ax + By + C
             *  So chi2 = SUM_i (Ax[i] + By[i] + C - u[i])^2
             *
             *  Must minimize wrt A, B & C simultaneously, i.e. solve the
             *  matrix equation
             *    M . X = Y
             *  where:
             *    X = (A,
             *         B,
             *         C)
             *
             *    Y = (SUM_i x_i u_i,
             *         SUM_i y_i u_i,
             *         SUM_i u_i)
             *  and
             *    M = ((SUM_i x_i ^ 2, SUM_i x_i y_i, SUM_i x_i),
             *         (SUM_i x_i y_i, SUM_i y_i ^ 2, SUM_i y_i),
             *         (SUM_i x_i,     SUM_i y_i,     N))
             *
             * So X = M_inv . Y
             * M is constant for a given virtual site so compute Minv once.
             *
             * In the case that M is (near-)singular, just take the average of
             * the connected points, i.e. replace M_inv with:
             *  (( 0, 0,   0),
             *   ( 0, 0,   0),
             *   ( 0, 0, 1/N)
             */

            // Compute y, by summing over neighbouring fluid sites.
            distribn_t sums[3] = { 0, 0, 0 };
            // {sumXU, sumYU, sumU}
            for (unsigned i = 0; i < vSite.neighbourGlobalIds.size(); ++i)
            {
              RSHV& hv = GetHV(latDat, hydroVarsCache, vSite.neighbourGlobalIds[i], t);
              LatticeSpeed uNorm = Dot(hv.u, iolet.GetNormal());
              sums[0] += hv.posIolet.x() * uNorm;
              sums[1] += hv.posIolet.y() * uNorm;
              sums[2] += uNorm;
            }

            // {A, B, C}
            util::Vector3D<distribn_t> coeffs = { 0, 0, 0 };

            // Do the matrix multiply to get the fitted coefficients.
            for (unsigned i = 0; i < 3; ++i)
              for (unsigned j = 0; j < 3; ++j)
                coeffs[i] += vSite.velocityMatrixInv[i][j] * sums[j];

            // Compute the magnitude of the velocity.
            LatticeSpeed ansNorm = Dot(coeffs, vSite.hv.posIolet);

            // multiply by the iolet normal and we're done!
            return iolet.GetNormal() * ansNorm;

          }

          RSHV& GetHV(const geometry::FieldData& latDat, RSHV::Map& hydroVarsCache,
                      const site_t globalIdx, const LatticeTimeStep t)
          {
            RSHV& ans = hydroVarsCache.at(globalIdx);
            /* Local sites have their entry in the cache set during collision
             * so they are guaranteed to be up to date. Neighbouring sites may
             * not be, but all the communication needed has been done. If the
             * calculation hasn't been done, do it and cache the results.
             */
            if (ans.t == t)
              return ans;

            auto neigh =
                latDat.GetNeighbouringData().GetSite(globalIdx);
            auto fOld = neigh.GetFOld<LatticeType>();
            LatticeType::CalculateDensityAndMomentum(fOld, ans.rho, ans.u);
            if (LatticeType::IsLatticeCompressible())
            {
              ans.u /= ans.rho;
            }
            ans.t = t;
            return ans;
          }
      };
    }

#endif /* HEMELB_LB_STREAMERS_VIRTUALSITEIOLET_H */
