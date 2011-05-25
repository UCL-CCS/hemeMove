#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "io/Writer.h"
#include "vis/ColPixel.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    class Screen
    {
      public:
        Screen();
        ~Screen();

        void
        AddPixel(const ColPixel* newPixel, const VisSettings* visSettings);
        void RenderLine(const float endPoint1[3],
                        const float endPoint2[3],
                        const VisSettings* visSettings);

        void Set(float maxX,
                 float maxY,
                 int pixelsX,
                 int pixelsY,
                 float rad,
                 const Viewpoint* viewpoint);

        void Resize(unsigned int pixelsX, unsigned int pixelsY);

        void Reset();

        /**
         * Does a transform from input array into output array. This function
         * will still work if the two arrays point to the same location in
         * memory. It only operates on the first two elements of input.
         *
         * @param input
         * @param output
         */
        template<typename T>
        void Transform(float* input, T output[2]) const
        {
          output[0] = (T) (ScaleX * (input[0] + MaxXValue));
          output[1] = (T) (ScaleY * (input[1] + MaxYValue));
        }

        const float* GetVtx() const;
        const float* GetUnitVectorProjectionX() const;
        const float* GetUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

        static const unsigned int COLOURED_PIXELS_MAX = 2048 * 2048;

        void
        CompositeImage(const VisSettings* visSettings);

        bool MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress);

        void WritePixelCount(io::Writer* writer);

        void WritePixels(const DomainStats* domainStats,
                         const VisSettings* visSettings,
                         io::Writer* writer);

        unsigned int GetPixelCount() const;

      private:
        template<bool xLimited> void RenderLineHelper(int x,
                                                      int y,
                                                      int incE,
                                                      int incNE,
                                                      int limit,
                                                      const VisSettings* visSettings);

        float ScaleX, ScaleY;
        float MaxXValue, MaxYValue;
        float vtx[3];

        // Projection of unit vectors along screen axes into normal space.
        float UnitVectorProjectionX[3];
        float UnitVectorProjectionY[3];

        unsigned int PixelsX, PixelsY;
        unsigned int PixelsMax;

        // Array of pixel ids.
        int* col_pixel_id;
        ColPixel localPixels[COLOURED_PIXELS_MAX];
        ColPixel compositingBuffer[COLOURED_PIXELS_MAX];
        // number of ColPixels.
        unsigned int col_pixels;
        unsigned int pixelCountInBuffer; // number received?
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
