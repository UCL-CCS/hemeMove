#ifndef HEMELB_VIS_CONTROL_H
#define HEMELB_VIS_CONTROL_H

#include <stack>

#include "geometry/LatticeData.h"

#include "lb/LbmParameters.h"
#include "lb/SimulationState.h"

#include "net/net.h"
#include "net/PhasedBroadcastIrregular.h"

#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

#include "vis/ColPixel.h"
#include "vis/GlyphDrawer.h"
#include "vis/StreaklineDrawer.h"
#include "vis/rayTracer/RayTracerEnhanced.h"

namespace hemelb
{
  namespace vis
  {
    /* the last three digits of the pixel identifier are used to
     indicate if the pixel is coloured via the ray tracing technique
     and/or a glyph and/or a particle/pathlet */

    /**
     * Class to control and use the effects of different visualisation methods.
     *
     * We use irregular phased broadcasting because we don't know in advance which iterations we'll
     * need to generate images on.
     *
     * The initial action is used to render the image on each core. 2 communications are required
     * between each pair of nodes so that the number of pixels can be communicated before the pixels
     * themselves. No overlap is possible between communications at different depths as the pixels
     * must be merged before they can be passed on. We don't need to pass info top-down, we only
     * pass image components upwards towards the top node.
     */
    class Control : public net::PhasedBroadcastIrregular<true, 2, 0, false, true>
    {
      public:
        Control(lb::StressTypes iStressType,
                net::Net* net,
                lb::SimulationState* simState,
                geometry::LatticeData* iLatDat);
        ~Control();

        void SetSomeParams(const float iBrightness,
                           const distribn_t iDensityThresholdMin,
                           const distribn_t iDensityThresholdMinMaxInv,
                           const distribn_t iVelocityThresholdMaxInv,
                           const distribn_t iStressThresholdMaxInv);

        void SetProjection(const int &pixels_x,
                           const int &pixels_y,
                           const float &ctr_x,
                           const float &ctr_y,
                           const float &ctr_z,
                           const float &longitude,
                           const float &latitude,
                           const float &zoom);

        bool MouseIsOverPixel(float* density, float* stress);

        void ProgressStreaklines(unsigned long time_step, unsigned long period);
        void Reset();

        void UpdateImageSize(int pixels_x, int pixels_y);
        void SetMouseParams(double iPhysicalPressure, double iPhysicalStress);
        void RegisterSite(site_t i, distribn_t density, distribn_t velocity, distribn_t stress);

        const ScreenPixels* GetResult(unsigned long startIteration);

        bool IsRendering() const;

        Viewpoint mViewpoint;
        DomainStats mDomainStats;
        VisSettings mVisSettings;

        double GetTimeSpent() const;

      protected:
        void InitialAction(unsigned long startIteration);
        void ProgressFromChildren(unsigned long startIteration, unsigned long splayNumber);
        void ProgressToParent(unsigned long startIteration, unsigned long splayNumber);
        void PostReceiveFromChildren(unsigned long startIteration, unsigned long splayNumber);
        void ClearOut(unsigned long startIteration);
        void InstantBroadcast(unsigned long startIteration);

      private:
        typedef net::PhasedBroadcastIrregular<true, 2, 0, false, true> base;
        typedef net::PhasedBroadcast<true, 2, 0, false, true> deepbase;
        typedef std::map<unsigned long, ScreenPixels*> mapType;

        // This is mainly constrained by the memory available per core.
        static const unsigned int SPREADFACTOR = 2;

        struct Vis
        {
            float half_dim[3];
            float system_size;
        };

        void initLayers();
        void Render();
        ScreenPixels* GetReceiveBuffer(unsigned int startIteration, unsigned int child);
        ScreenPixels* GetPixFromBuffer();

        // Because of the 2-splay, we need to have two sets of receive buffers, so that comms
        // on consecutive iterations don't overwrite one another.
        ScreenPixels recvBuffers[2][SPREADFACTOR];

        std::stack<ScreenPixels*> pixelsBuffer;
        std::map<unsigned long, ScreenPixels*> resultsByStartIt;
        geometry::LatticeData* mLatDat;
        Screen mScreen;
        Vis* vis;
        raytracer::RayTracerEnhanced *myRayTracer;
        GlyphDrawer *myGlypher;
        StreaklineDrawer *myStreaker;

        double timeSpent;
    };
  }
}

#endif // HEMELB_VIS_CONTROL_H
