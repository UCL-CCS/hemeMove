
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_BUILDSYSTEMINTERFACE_H
#define HEMELB_LB_BUILDSYSTEMINTERFACE_H

#include "lb/kernels/Kernels.h"
#include "lb/kernels/momentBasis/MomentBases.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "lb/collisions/Collisions.h"
#include "lb/streamers/Streamers.h"

namespace hemelb
{
  namespace lb
  {
    /**
     * The names of the following classes must correspond to options given for the CMake
     * HEMELB_KERNEL parameter.
     */
    template<class Lattice>
    class LBGK
    {
      public:
        typedef kernels::LBGK<Lattice> Type;
    };

    template<class Lattice>
    class ADVECTIONDIFFUSIONLBGK
    {
      public:
        typedef kernels::AdvectionDiffusionLBGK<Lattice> Type;
    };

    /**
     * The entropic implementation by Ansumali et al.
     */
    template<class Lattice>
    class EntropicAnsumali
    {
      public:
        typedef kernels::EntropicAnsumali<Lattice> Type;
    };

    /**
     * The entropic implementation by Chikatamarla et al.
     */
    template<class Lattice>
    class EntropicChik
    {
      public:
        typedef kernels::EntropicChik<Lattice> Type;
    };

    /**
     * MRT currently we only have DHumieres implementation, on D3Q15 and D3Q19 lattices.
     */
    template<class Lattice>
    class MRT
    {
    };

    template<class Lattice>
    class ADVECTIONDIFFUSIONMRT
    {
    };

    template<>
    class MRT<lattices::D3Q15>
    {
      public:
        typedef kernels::MRT<kernels::momentBasis::DHumieresD3Q15MRTBasis> Type;
    };

    template<>
    class ADVECTIONDIFFUSIONMRT<lattices::D3Q15>
    {
      public:
        typedef kernels::AdvectionDiffusionMRT<kernels::momentBasis::DHumieresD3Q15MRTBasis> Type;
    };

    template<>
    class MRT<lattices::D3Q19>
    {
      public:
        typedef kernels::MRT<kernels::momentBasis::DHumieresD3Q19MRTBasis> Type;
    };

    template<>
    class ADVECTIONDIFFUSIONMRT<lattices::D3Q19>
    {
      public:
        typedef kernels::AdvectionDiffusionMRT<kernels::momentBasis::DHumieresD3Q19MRTBasis> Type;
    };

    template<class Lattice>
    class TRT
    {
      public:
        typedef kernels::TRT<Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with Carreau-Yasuda rheology model.
     */
    template<class Lattice>
    class NNCY
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::CarreauYasudaRheologyModelHumanFit, Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with Carreau-Yasuda rheology model fitted to experimental data on murine blood viscosity.
     */
    template<class Lattice>
    class NNCYMOUSE
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::CarreauYasudaRheologyModelMouseFit, Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with Casson rheology model.
     */
    template<class Lattice>
    class NNC
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::CassonRheologyModel, Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with truncated power law rheology model.
     */
    template<class Lattice>
    class NNTPL
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::TruncatedPowerLawRheologyModel, Lattice> Type;
    };

    /**
     * The following classes have names corresponding to the options given in the build system for
     * HEMELB_WALL_BOUNDARY
     */
    /**
     * The Bouzidi-Firdaous-Lallemand interpolation-based boundary condition.
     */
    template<class Collision>
    class BFL
    {
      public:
        typedef typename streamers::BouzidiFirdaousLallemand<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLBFLIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallBFLIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLBFLIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallBFLIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLSBBIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallSBBIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLSBBIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallSBBIoletNeumann<Collision>::Type Type;
    };
    /**
     * The Guo Zheng and Shi mode-extrapolation boundary condition.
     */
    template<class Collision>
    class GZS
    {
      public:
        typedef typename streamers::GuoZhengShi<Collision>::Type Type;
    };
    /**
     * The simple bounce back boundary condition.
     */
    template<class Collision>
    class SIMPLEBOUNCEBACK
    {
      public:
        typedef typename streamers::SimpleBounceBack<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLBFLIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallBFLIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLBFLIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallBFLIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLSBBIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallSBBIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLSBBIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallSBBIoletNeumann<Collision>::Type Type;
    };
    /**
     * The Junk & Yang 2005 boundary condition.
     */
    template<class Collision>
    class JUNKYANG
    {
      public:
        typedef typename streamers::JunkYang<Collision>::Type Type;
    };

    /**
     * The following classes have names corresponding to the options given in the build system for
     * HEMELB_INLET_BOUNDARY / HEMELB_OUTLET_BOUNDARY
     */

    /**
     * Our zeroth-order phantom site BC for iolets
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREIOLET
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIolet<Collision>::Type Type;
    };
    /**
     * The inlet/outlet condition based on Ladd's modified bounce-back on
     * links.
     */
    template<class Collision>
    struct LADDIOLET
    {
        typedef typename streamers::LaddIolet<Collision>::Type Type;
    };

    /**
     * The following classes have names corresponding to the options given in the build system for
     * HEMELB_WALL_INLET_BOUNDARY / HEMELB_WALL_OUTLET_BOUNDARY
     */
    /**
     * Nash in/outlet + SBB
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSURESBB
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIoletSBB<Collision>::Type Type;
    };

    /**
     * Ladd in/outlet + SBB
     */
    template<class Collision>
    struct LADDIOLETSBB
    {
        typedef typename streamers::LaddIoletSBB<Collision>::Type Type;
    };

    /**
     * Nash in/outlet + BFL
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREBFL
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIoletBFL<Collision>::Type Type;
    };

    /**
     * Ladd in/outlet + BFL
     */
    template<class Collision>
    struct LADDIOLETBFL
    {
        typedef typename streamers::LaddIoletBFL<Collision>::Type Type;
    };
    /**
     * Nash in/outlet + GZS
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREGZS
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIoletGZS<Collision>::Type Type;
    };

    /**
     * Ladd in/outlet + GZS
     */
    template<class Collision>
    struct LADDIOLETGZS
    {
        typedef typename streamers::LaddIoletGZS<Collision>::Type Type;
    };

    /**
     * Nash/Krueger in/outlet + SBB
     */
    template<class Collision>
    struct VIRTUALSITEIOLETSBB
    {
        typedef typename streamers::VirtualSiteIolet<Collision> Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONOUTFLOWIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionOutflowIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONINLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionInletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLINLETDIRICHLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallInletDirichletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLINLETDIRICHLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallInletDirichletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLINLETDIRICHLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallInletDirichletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLINLETDIRICHLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallInletDirichletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONOUTFLOWBOUNCEBACKIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionOutflowBounceBackIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWBOUNCEBACKIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowBounceBackIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWBOUNCEBACKIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowBounceBackIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWBOUNCEBACKIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowBounceBackIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWBOUNCEBACKIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowBounceBackIoletNeumann<Collision>::Type Type;
    };
  }
}

#endif /* HEMELB_LB_BUILDSYSTEMINTERFACE_H */
