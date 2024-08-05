// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LB_HPP
#define HEMELB_LB_LB_HPP

#include "lb/lb.h"
#include <concepts>
#include "lb/InitialCondition.h"
#include "lb/InitialCondition.hpp"

namespace hemelb::lb
{

    namespace detail {
        // Helper for below using the standard trick of deducing the
        // indices from an index_sequence.
        template<typename F, typename... Ts, std::size_t... Is>
        void tuple_enumerated_foreach_impl(F &&f, std::tuple<Ts...>& t, std::index_sequence<Is...> is) {
            (std::invoke(f, Is, std::get<Is>(t)), ...);
        }
    }

    // Invoke the function object with arguments (i, std::get<i>(tuple)),
    // where i is a std::size_t and the tuple is the other argument.
    //
    // Since element types of a tuple can be heterogeneous, the function
    // object needs to be similar to a variant's visitor and have an
    // overloaded or generic operator().
    //
    // TODO: consider a const overload
    template <typename F, typename... Ts>
    void tuple_enumerated_foreach(std::tuple<Ts...>& t, F&& func) {
        using I = std::make_index_sequence<sizeof...(Ts)>;
        detail::tuple_enumerated_foreach_impl<F, Ts...>(std::forward<F>(func), t, I{});
    }

    template<class TRAITS>
    void LBM<TRAITS>::InitCollisions()
    {
        // Ensure the boundary objects have all info necessary.
        PrepareBoundaryObjects();

        // TODO Note that the convergence checking is not yet implemented in the
        // new boundary condition hierarchy system.
        // It'd be nice to do this with something like
        // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

        auto& dom = mLatDat->GetDomain();
        auto initParams = InitParams();
        initParams.latDat = &dom;
        initParams.lbmParams = &mParams;
        initParams.neighbouringDataManager = neighbouringDataManager;
        initParams.siteRanges.resize(2);

        auto get_boundary_vals = [this](std::size_t collision_id) -> BoundaryValues* {
            if (collision_id < 2)
                return nullptr;
            if (collision_id % 2)
                return mOutletValues;
            else
                return mInletValues;
        };

        tuple_enumerated_foreach(
                mCollisions,
                [&] <typename COLL> (std::size_t i, std::unique_ptr<COLL>& collision) {
                    initParams.siteRanges[0] = dom.GetMidDomainSiteRange(i);
                    initParams.siteRanges[1] = dom.GetDomainEdgeSiteRange(i);
                    initParams.boundaryObject = get_boundary_vals(i);
                    collision = std::make_unique<COLL>(initParams);
                }
        );
    }



    template<class TRAITS>
    void LBM<TRAITS>::SetInitialConditions(lb::InitialCondition const& icond, const net::IOCommunicator& ioComms)
    {
      icond.SetFs<LatticeType>(mLatDat, ioComms);
      icond.SetTime(mState);
    }

    template<class TRAITS>
    void LBM<TRAITS>::PreSend()
    {
        timings.lb().Start();
        timings.lb_calc().Start();

        // In the PreSend phase, we do LB on all the sites that need to have results sent to
        // neighbouring ranks ('domainEdge' sites). In site id terms, this means we start at the
        // end of the sites whose neighbours all lie on this rank ('midDomain'), then progress
        // through the sites of each type in turn.
        log::Logger::Log<log::Debug, log::OnePerCore>("LBM - PreSend - StreamAndCollide");
        tuple_enumerated_foreach(
                mCollisions,
                [this] (std::size_t i, auto& collision) {
                    auto [begin, end] = mLatDat->GetDomain().GetDomainEdgeSiteRange(i);
                    StreamAndCollide(*collision, begin, end);
                }
        );

        timings.lb_calc().Stop();
        timings.lb().Stop();
    }

    template<class TRAITS>
    void LBM<TRAITS>::PreReceive()
    {
        timings.lb().Start();
        timings.lb_calc().Start();

        // In the PreReceive phase, we perform LB for all the sites whose neighbours lie on this
        // rank ('midDomain' rather than 'domainEdge' sites). Ideally this phase is the longest bit (maximising time for the asynchronous sends
        // and receives to complete).
        //
        // In site id terms, this means starting at the first site and progressing through the
        // midDomain sites, one type at a time.
        log::Logger::Log<log::Debug, log::OnePerCore>("LBM - PreReceive - StreamAndCollide");
        tuple_enumerated_foreach(
                mCollisions,
                [this] (std::size_t i, auto& collision) {
                    auto [begin, end] = mLatDat->GetDomain().GetMidDomainSiteRange(i);
                    StreamAndCollide(*collision, begin, end);
                }
        );
        timings.lb_calc().Stop();
        timings.lb().Stop();
    }

    template<class TRAITS>
    void LBM<TRAITS>::PostReceive()
    {
        timings.lb().Start();

        // Copy the distribution functions received from the neighbouring
        // processors into the destination buffer "f_new".
        // This is done here, after receiving the sent distributions from neighbours.
        mLatDat->CopyReceived();

        // Do any cleanup steps necessary on boundary nodes
        timings.lb_calc().Start();

        log::Logger::Log<log::Debug, log::OnePerCore>("LBM - PostReceive - StreamAndCollide");
        tuple_enumerated_foreach(
                mCollisions,
                [this] (std::size_t i, auto& collision) {
                    auto [begin, end] = mLatDat->GetDomain().GetDomainEdgeSiteRange(i);
                    PostStep(*collision, begin, end);
                }
        );
        tuple_enumerated_foreach(
                mCollisions,
                [this] (std::size_t i, auto& collision) {
                    auto [begin, end] = mLatDat->GetDomain().GetMidDomainSiteRange(i);
                    PostStep(*collision, begin, end);
                }
        );

        timings.lb_calc().Stop();
        timings.lb().Stop();
    }
}

#endif /* HEMELB_LB_LB_HPP */
