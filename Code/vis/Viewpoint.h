#ifndef HEMELB_VIS_VIEWPOINT_H
#define HEMELB_VIS_VIEWPOINT_H

namespace hemelb
{
  namespace vis
  {
    class Viewpoint
    {
      public:
        void RotateToViewpoint(float iXIn,
                               float iYIn,
                               float iZIn,
                               float *oXOut,
                               float *oYOut,
                               float *oZOut);

        void Project(const float p1[], float p2[]);

        float x[3];
        float SinYRotation, CosYRotation;
        float SinXRotation, CosXRotation;
        float dist;

      private:
        void Rotate(float sinX,
                    float cosX,
                    float sinY,
                    float cosY,
                    float xIn,
                    float yIn,
                    float zIn,
                    float* xOut,
                    float* yOut,
                    float* zOut);

    };
  }
}

#endif /* HEMELB_VIS_VIEWPOINT_H */
