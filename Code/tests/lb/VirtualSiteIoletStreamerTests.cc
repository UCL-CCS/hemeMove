// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

#include "lb/Kernels.h"
#include "lb/Streamers.h"
#include "lb/lattices/D3Q15.h"
#include "lb/iolets/InOutLets.h"
#include "lb/streamers/VirtualSiteIolet.h"
#include "geometry/SiteData.h"
#include "util/numerical.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb::tests
{
    constexpr distribn_t allowedError = 1e-10;
    static Approx apprx(double x) {
        return Approx(x).margin(allowedError);
    }

    /**
     * Helper class that exposes implementation details to the tests.
     */
    using lb::VSExtra;
    class LocalVSExtra : public VSExtra<lb::D3Q15>
    {
    public:
        //typedef hemelb::lb::streamers::VSExtra::UnitVec UnitVec;
        LocalVSExtra(lb::InOutLet& iolet) :
                VSExtra<lb::D3Q15> (iolet)
        {
        }
        [[nodiscard]] const UnitVec& GetE1() const
        {
            return e1;
        }
        [[nodiscard]] const UnitVec& GetE2() const
        {
            return e2;
        }
        [[nodiscard]] const UnitVec& GetE3() const
        {
            return n;
        }
    };

    static LatticeVelocity GetVelocity(LatticeVector pos)
    {
        LatticeVelocity c(0.01, 0.01, 0.);
        LatticeVelocity u(0.);
        u.z() = Dot(c, pos);
        return u;
    }
    static LatticeDensity GetDensity(LatticeVector pos)
    {
        util::Vector3D<LatticeDensity> gradRho(0., 0., -1e-2);
        return 1.045 + Dot(gradRho, pos);
    }

    using lb::InOutLetCosine;
    static InOutLetCosine* GetIolet(lb::BoundaryValues& iolets)
    {
        auto ans = dynamic_cast<lb::InOutLetCosine*> (iolets.GetLocalIolet(0));
        REQUIRE(ans != nullptr);
        return ans;
    }

    /**
     * VirtualSiteIoletStreamerTests:
     *
     * This class tests the streamer implementations. We assume the
     * collision operators are correct (as they're tested elsewhere),
     * then compare the post-streamed values with the values we expect
     * to have been streamed there.
     */
    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "VirtualSiteIoletStreamerTests") {

        using Lattice = lb::D3Q15;
        using Kernel = lb::LBGK<Lattice>;
        using Collision = lb::Normal<Kernel>;
        using VirtualSite = lb::VirtualSite<Lattice>;
        using InOutLetCosine = lb::InOutLetCosine;
        using RSHV = lb::RSHV;

        auto Info = [&]() -> const lb::LatticeInfo& {
            return Lattice::GetLatticeInfo();
        };

        auto CheckAllHVUpdated = [&](lb::BoundaryValues& iolets, LatticeTimeStep expectedT) {
            auto * extra = dynamic_cast<VSExtra<Lattice>*> (iolets.GetLocalIolet(0)->GetExtraData());
            for (auto& [siteGlobalIdx, hv] : extra->hydroVarsCache) {
                LatticeVector sitePos;
                dom->GetGlobalCoordsFromGlobalNoncontiguousSiteId(siteGlobalIdx, sitePos);
                REQUIRE(expectedT == hv.t);
                REQUIRE(apprx(GetDensity(sitePos)) == hv.rho);
                LatticeVelocity u = GetVelocity(sitePos);
                for (unsigned i = 0; i < 3; ++i)
                    REQUIRE(apprx(u[i]) == hv.u[i]);
            }
        };

        auto CheckPointInCache = [&](LocalVSExtra& extra, LatticeVector expectedPt,
                                     LatticePosition expectedIoletPos)
        {
            site_t expectedGlobalIdx =
                    dom->GetGlobalNoncontiguousSiteIdFromGlobalCoords(expectedPt);
            RSHV::Map::iterator hvPtr = extra.hydroVarsCache.find(expectedGlobalIdx);

            REQUIRE(hvPtr != extra.hydroVarsCache.end());
            for (unsigned i = 0; i < 3; ++i)
                REQUIRE(apprx(expectedIoletPos[i]) == hvPtr->second.posIolet[i]);
        };

        auto InitialiseGradientHydroVars = [&](){
            for (unsigned i = 1; i <= 4; ++i) {
                for (unsigned j = 1; j <= 4; ++j) {
                    for (unsigned k = 1; k <= 4; ++k) {
                        LatticeVector pos(i, j, k);
                        site_t siteIdx = dom->GetContiguousSiteId(pos);
                        //geometry::Site < geometry::domain_type > site = latDat->GetSite(siteIdx);
                        auto fOld = Lattice::mut_span{latDat->GetFNew(siteIdx * Lattice::NUMVECTORS), Lattice::NUMVECTORS};
                        LatticeDensity rho = GetDensity(pos);
                        LatticeVelocity u = GetVelocity(pos);
                        u *= rho;
                        Lattice::CalculateFeq(rho, u, fOld);
                    }
                }
            }
            latDat->SwapOldAndNew();
        };

        auto propertyCache = lb::MacroscopicPropertyCache(*simState, *dom);

        auto inletBoundary = BuildIolets(geometry::INLET_TYPE);
        InOutLetCosine* inlet = GetIolet(inletBoundary);
        // We have to make the outlet sane and consistent with the geometry now.
        inlet->SetNormal(util::Vector3D<Dimensionless>(0, 0, 1));
        PhysicalPosition inletCentre(2.5, 2.5, 0.5);
        inletCentre *= simConfig.GetVoxelSize();
        inlet->SetPosition(unitConverter->ConvertPositionToLatticeUnits(inletCentre));
        // Want to set the density gradient to be 0.01 in lattice units,
        // starting at 1.0 at the outlet.
        inlet->SetPressureAmp(0.);
        inlet->SetPressureMean(unitConverter->ConvertPressureToPhysicalUnits(1.04 * Cs2));

        // Same for the cut distances of outlet sites.
        for (unsigned i = 1; i <= 4; ++i) {
            for (unsigned j = 1; j <= 4; ++j) {
                site_t siteId = dom->GetContiguousSiteId(LatticeVector(i, j, 1));
                geometry::Site < geometry::Domain > site = latDat->GetSite(siteId);
                for (Direction p = 0; p < Info().GetNumVectors(); ++p) {
                    if (Info().GetVector(p).z() < 0.) {
                        // Sanity check
                        REQUIRE(site.HasIolet(p));
                        // Set the cut distance to half way
                        dom->SetBoundaryDistance(siteId, p, 0.5);
                    }
                }
            }
        }

        auto outletBoundary = BuildIolets(geometry::OUTLET_TYPE);

        InOutLetCosine* outlet = GetIolet(outletBoundary);
        // We have to make the outlet sane and consistent with the geometry now.
        outlet->SetNormal(util::Vector3D<Dimensionless>(0, 0, -1));
        PhysicalPosition outletCentre(2.5, 2.5, 4.5);
        outletCentre *= simConfig.GetVoxelSize();
        outlet->SetPosition(unitConverter->ConvertPositionToLatticeUnits(outletCentre));
        outlet->SetPressureAmp(0.);
        outlet->SetPressureMean(unitConverter->ConvertPressureToPhysicalUnits(1.0 * Cs2));

        // Same for the cut distances of outlet sites.
        for (unsigned i = 1; i <= 4; ++i) {
            for (unsigned j = 1; j <= 4; ++j) {
                site_t siteId = dom->GetContiguousSiteId(LatticeVector(i, j, 4));
                geometry::Site < geometry::Domain > site = latDat->GetSite(siteId);
                for (Direction p = 0; p < Info().GetNumVectors(); ++p) {
                    if (Info().GetVector(p).z() > 0.) {
                        // Sanity check
                        REQUIRE(site.HasIolet(p));
                        // Set the cut distance to half way
                        dom->SetBoundaryDistance(siteId, p, 0.5);
                    }
                }
            }
        }

        SECTION("TestVirtualSiteMatrixInverse") {
            distribn_t inv[3][3];
            distribn_t det;

            // Try the identity matrix first
            distribn_t eye[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            det = VirtualSite::Matrix3DInverse(eye, inv);
            REQUIRE(apprx(1.) == det);
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    REQUIRE(apprx(eye[i][j]) == inv[i][j]);

            // arbitrary matrix with nice ish values
            distribn_t a[3][3] = { { 1, 2, 3 }, { 0, 4, 5 }, { 1, 0, 6 } };
            distribn_t aInv22[3][3] = { { 24, -12, -2 }, { 5, 3, -5 }, { -4, 2, 4 } };

            det = VirtualSite::Matrix3DInverse(a, inv);
            REQUIRE(apprx(22.) == det);
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    REQUIRE(apprx(aInv22[i][j] / 22) == inv[i][j]);

        }

        SECTION("TestVirtualSiteConstruction") {
            InOutLetCosine& outlet = *GetIolet(outletBoundary);
            LocalVSExtra extra(outlet);
            // Expected basis for iolet coords
            util::Vector3D<Dimensionless> e1(0, 1, 0), e2(1, 0, 0), e3(0, 0, -1);
            // I.e. flip x & y, negate z.
            REQUIRE(ApproxV(e1).Margin(allowedError) == extra.GetE1());
            REQUIRE(ApproxV(e2).Margin(allowedError) == extra.GetE2());
            REQUIRE(ApproxV(e3).Margin(allowedError) == extra.GetE3());

            // This point should be just outside the 4 cube's outlet
            LatticeVector vLoc(2, 2, 5);
            // Create it
            VirtualSite vSite(initParams, extra, vLoc);

            // The iolet coords
            // position relative to outlet (world coords) = (-0.5, -0.5, 0.5)
            // Expect in iolet coords this to be (-0.5, -0.5, -0.5)
            REQUIRE(ApproxVector<LatticePosition>(-0.5, -0.5, -0.5).Margin(allowedError) == vSite.hv.posIolet);

            // For D3Q15, there should be 5 cut links.
            REQUIRE(size_t(5) == vSite.neighbourGlobalIds.size());
            // each with q = 0.5
            REQUIRE(apprx(1.25) == vSite.sumQiSq);
            // For each site index in the neighbourhood
            for (Direction p = 0; p < Info().GetNumVectors(); ++p) {
//              if (Info().GetVector(p).z < 0.)
//              {
//                CPPUNIT_ASSERT(vSite.streamingIndices[p] < 4 * 4 * 4 * 15);
//              }
//              else
//              {
//                CPPUNIT_ASSERT_EQUAL(site_t(4 * 4 * 4 * 15), vSite.streamingIndices[p]);
//              }
            }

            // Now check that appropriate entries have been added to the hydroCache
            CheckPointInCache(extra, LatticeVector(1, 1, 4), LatticePosition(-1.5, -1.5, 0.5));
            CheckPointInCache(extra, LatticeVector(1, 3, 4), LatticePosition(0.5, -1.5, 0.5));
            CheckPointInCache(extra, LatticeVector(3, 1, 4), LatticePosition(-1.5, 0.5, 0.5));
            CheckPointInCache(extra, LatticeVector(3, 3, 4), LatticePosition(0.5, 0.5, 0.5));
            CheckPointInCache(extra, LatticeVector(2, 2, 4), LatticePosition(-0.5, -0.5, 0.5));

            // The velocity matrix = sum over i of:
            // x_i^2   x_i y_i x_i
            // x_i y_i y_i^2   y_i
            // x_i     y_i     1
            // distribn_t velMat[3][3] = { { 5.25, 1.25, -2.5 },
            //                             { 1.25, 5.25, -2.5 },
            //                             { -2.5, -2.5,  5.0 } };
            // It's inverse is
            distribn_t velMatInv[3][3] = { { 0.25, 0., 0.125 },
                                           { 0., 0.25, 0.125 },
                                           { 0.125, 0.125, 0.325 } };
            for (unsigned i = 0; i < 3; ++i) {
                for (unsigned j = 0; j < 3; ++j) {
                    REQUIRE(apprx(velMatInv[i][j]) == vSite.velocityMatrixInv[i][j]);
                }
            }

        }

        SECTION("TestStreamerInitialisation") {
            initParams.boundaryObject = &outletBoundary;
            // Set up the ranges to cover Mid 3 (pure outlet) and Mid 5 (outlet/wall)
            initParams.siteRanges = {dom->GetMidDomainSiteRange(3), dom->GetMidDomainSiteRange(5)};
            lb::VirtualSiteIolet<Collision> outletStreamer(initParams);

            // All the sites at the outlet plane (x, y, 3) should be in the cache.
            InOutLetCosine& outlet = *GetIolet(outletBoundary);
            auto* extra = dynamic_cast<VSExtra<Lattice>*> (outlet.GetExtraData());
            REQUIRE(extra != nullptr);

            for (unsigned i = 1; i <= 4; ++i) {
                for (unsigned j = 1; j <= 4; ++j) {
                    LatticeVector pos(i, j, 4);
                    site_t globalIdx = dom->GetGlobalNoncontiguousSiteIdFromGlobalCoords(pos);
                    //                site_t localIdx = latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(globalIdx);
                    //                geometry::Site < geometry::domain_type > site = latDat->GetSite(localIdx);

                    RSHV::Map::iterator hvPtr = extra->hydroVarsCache.find(globalIdx);
                    REQUIRE(hvPtr != extra->hydroVarsCache.end());
                }
            }

            // And the reverse is true: every cache entry should be a site at the outlet plane
            for (auto& hvPtr : extra->hydroVarsCache) {
                site_t globalIdx = hvPtr.first;
                LatticeVector pos;
                dom->GetGlobalCoordsFromGlobalNoncontiguousSiteId(globalIdx, pos);
                REQUIRE(util::IsInRange<LatticeCoordinate>(pos.x(), 1, 4));
                REQUIRE(util::IsInRange<LatticeCoordinate>(pos.y(), 1, 4));
                REQUIRE(LatticeCoordinate(4) == pos.z());
            }
        }

        SECTION("TestStep") {
            initParams.boundaryObject = &inletBoundary;
            lb::VirtualSiteIolet<Collision> inletStreamer(initParams);
            initParams.boundaryObject = &outletBoundary;
            lb::VirtualSiteIolet<Collision> outletStreamer(initParams);

            InitialiseGradientHydroVars();

            // Stream and collide
            inletStreamer.StreamAndCollide(dom->GetMidDomainSiteRange(2).first,
                                           dom->GetMidDomainSiteRange(2).second,
                                           &lbmParams,
                                           *latDat ,
                                           propertyCache);

            outletStreamer.StreamAndCollide(dom->GetMidDomainSiteRange(3).first,
                                            dom->GetMidDomainSiteRange(3).second,
                                            &lbmParams,
                                            *latDat,
                                            propertyCache);

            inletStreamer.StreamAndCollide(dom->GetMidDomainSiteRange(4).first,
                                           dom->GetMidDomainSiteRange(4).second,
                                           &lbmParams,
                                           *latDat,
                                           propertyCache);

            outletStreamer.StreamAndCollide(dom->GetMidDomainSiteRange(5).first,
                                            dom->GetMidDomainSiteRange(5).second,
                                            &lbmParams,
                                            *latDat,
                                            propertyCache);

            // Now every entry in the RSHV cache should have been updated
            CheckAllHVUpdated(inletBoundary, 0);
            CheckAllHVUpdated(outletBoundary, 0);

            // Stream and collide
            inletStreamer.PostStep(dom->GetMidDomainSiteRange(2).first,
                                   dom->GetMidDomainSiteRange(2).second,
                                   &lbmParams,
                                   latDat.get(),
                                   propertyCache);

            outletStreamer.PostStep(dom->GetMidDomainSiteRange(3).first,
                                    dom->GetMidDomainSiteRange(3).second,
                                      &lbmParams,
                                      latDat.get(),
                                      propertyCache);

            inletStreamer.PostStep(dom->GetMidDomainSiteRange(4).first,
                                   dom->GetMidDomainSiteRange(4).second,
                                     &lbmParams,
                                     latDat.get(),
                                     propertyCache);

            outletStreamer.PostStep(dom->GetMidDomainSiteRange(5).first,
                                    dom->GetMidDomainSiteRange(5).second,
                                      &lbmParams,
                                      latDat.get(),
                                      propertyCache);

            // Check that all the vsites have sensible hydro values
            InOutLetCosine* inlet = GetIolet(inletBoundary);
            auto * inExtra = dynamic_cast<VSExtra<Lattice>*> (inlet->GetExtraData());

            for (auto& [vSiteGlobalIdx, vSite] : inExtra->vSites) {
                REQUIRE(LatticeTimeStep(1) == vSite.hv.t);
                REQUIRE(apprx(LatticeDensity(1.045)) == vSite.hv.rho);

                LatticeVector pos;
                dom->GetGlobalCoordsFromGlobalNoncontiguousSiteId(vSiteGlobalIdx, pos);

                if (vSite.neighbourGlobalIds.size() > 3)
                    REQUIRE(apprx(GetVelocity(pos).z()) ==  vSite.hv.u.z());
            }

            InOutLetCosine* outlet = GetIolet(outletBoundary);
            auto outExtra = dynamic_cast<VSExtra<Lattice>*> (outlet->GetExtraData());

            for (auto& [vSiteGlobalIdx, vSite] : outExtra->vSites) {
                REQUIRE(LatticeTimeStep(1) == vSite.hv.t);
                REQUIRE(apprx(LatticeDensity(0.995)) == vSite.hv.rho);

                LatticeVector pos;
                dom->GetGlobalCoordsFromGlobalNoncontiguousSiteId(vSiteGlobalIdx, pos);

                if (vSite.neighbourGlobalIds.size() > 3)
                    REQUIRE(apprx(GetVelocity(pos).z()) == vSite.hv.u.z());

            }
        }
    }
}
