#ifndef HEMELB_VIS_VIEWPOINT_H
#define HEMELB_VIS_VIEWPOINT_H

#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    /**
     * Holds the position of the viewpoint or camera location and performs 
     * projection of points in world coordinates into 3D screen co-ordinates
     */
    class Viewpoint
    {
        //See http://en.wikipedia.org/wiki/3D_projection#Perspective_projection
      public:
        Viewpoint();

        /**
         * Chanes a location in camera coordinates into world coordinates
         * Note that both co-ordinate systems share the same (0,0,0) point
         * and the camera is at (0,0,radius) in camera coordinates
         *
         */
        Vector3D<distribn_t>
        RotateCameraCoordinatesToWorldCoordinates(const Vector3D<distribn_t>& iVector) const;

        /**
         * Does the reverse of the above
         *
         */
        Vector3D<distribn_t> RotateWorldToCameraCoordinates(const Vector3D<distribn_t>& iVector) const;

        /**
         * Projects a location in world-cordinates onto the infinite screen,
         * by translating and rotating such as to give coordintes relative
         * to the camera, and then applying a perspective projection
         */
        Vector3D<distribn_t> Project(const Vector3D<distribn_t>& p1) const;

        /**
         * Sets the position of the Camera or Viewpoint
         *
         * In Word co-ordinates, the camera is located at radius from
         * the local centre, and pointing at an angle indicated by the
         * latitude and longitude.
         *
         * @param iLongitude - longitude in radians.
         * @param iLatitude - latitude in radians.
         * @param iLocalCentre - where the camera should point in world-coordinates
         * @param iRadius - the distance of the camera from this local centre
         * @param iDistanceFromEyeToScreen - the distance of the inifite screen
         * from the viewer. Allows for zoom
         */
        void SetViewpointPosition(float iLongitude,
                                  float iLatitude,
                                  const Vector3D<distribn_t>& iLocalCentre,
                                  float iRadius,
                                  float iDistanceFromEyeToScreen);

        const Vector3D<distribn_t>& GetViewpointCentreLocation() const;

      private:
        /**
         * Performs a vector rotation using stored
         * Sin and Cosine Values
         */
        static Vector3D<distribn_t> Rotate(float iSinThetaX,
                                           float iCosThetaX,
                                           float iSinThetaY,
                                           float iCosThetaY,
                                           const Vector3D<distribn_t>& iVector);

        /*
         * Reverses a vector rotation of the above
         */
        static Vector3D<distribn_t> UnRotate(float iSinThetaX,
                                             float iCosThetaX,
                                             float iThetaY,
                                             float iCosThetaY,
                                             const Vector3D<distribn_t>& iVector);

        float mDistanceFromEyeToScreen;

        float mSinLongitude;
        float mCosLongitude;
        float mSinLatitude;
        float mCosLatitude;

        // Stores the viewpoint Location in world co-ordinate
        Vector3D<distribn_t> mViewpointLocation;
    };
  }
}

#endif /* HEMELB_VIS_VIEWPOINT_H */
