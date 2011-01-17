#include <vector>
#include <math.h>

#include "util/utilityFunctions.h"
#include "vis/Control.h"
#include "vis/RayTracer.h"
#include "vis/GlyphDrawer.h"

#include "io/XdrFileWriter.h"

namespace hemelb
{
  namespace vis
  {
    // TODO ye gods.
    // make a global controller
    Control *controller;

    Control::Control(lb::StressTypes iStressType,
                     lb::GlobalLatticeData &iGlobLatDat)
    {
      mStressType = iStressType;

      this->vis = new Vis;

      //sites_x etc are globals declared in net.h
      vis->half_dim[0] = 0.5F * float(iGlobLatDat.GetXSiteCount());
      vis->half_dim[1] = 0.5F * float(iGlobLatDat.GetYSiteCount());
      vis->half_dim[2] = 0.5F * float(iGlobLatDat.GetZSiteCount());

      vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1],
                                                             vis->half_dim[2]));
      col_pixels_max = COLOURED_PIXELS_MAX;

      col_pixel_recv[0] = new ColPixel[col_pixels_max];
      col_pixel_recv[1] = new ColPixel[col_pixels_max];

      pixels_max = COLOURED_PIXELS_MAX;
      col_pixel_id = new int[pixels_max];

      for (int i = 0; i < COLOURED_PIXELS_MAX; i++)
      {
        col_pixel_id[i] = -1;
      }

    }

    void Control::initLayers(topology::NetworkTopology * iNetworkTopology,
                             lb::GlobalLatticeData &iGlobLatDat,
                             lb::LocalLatticeData &iLocalLatDat)
    {
      myRayTracer
          = new RayTracer(iNetworkTopology, &iLocalLatDat, &iGlobLatDat);
      myGlypher = new GlyphDrawer(&iGlobLatDat, &iLocalLatDat);

#ifndef NO_STREAKLINES
      myStreaker = new StreaklineDrawer(iNetworkTopology, iLocalLatDat,
                                        iGlobLatDat);
#endif
      // Note that rtInit does stuff to this->ctr_x (because this has
      // to be global)
      this->ctr_x -= vis->half_dim[0];
      this->ctr_y -= vis->half_dim[1];
      this->ctr_z -= vis->half_dim[2];
    }

    void Control::RotateAboutXThenY(const float &iSinThetaX,
                                    const float &iCosThetaX,
                                    const float &iSinThetaY,
                                    const float &iCosThetaY,
                                    const float &iXIn,
                                    const float &iYIn,
                                    const float &iZIn,
                                    float &oXOut,
                                    float &oYOut,
                                    float &oZOut)
    {
      // A rotation of iThetaX clockwise about the x-axis
      // Followed by a rotation of iThetaY anticlockwise about the y-axis.
      // In matrices:
      //       (cos(iThetaY)  0 sin(iThetaY)) (1 0            0              )
      // Out = (0             1 0           ) (0 cos(-iThetaX) -sin(-iThetaX)) In
      //       (-sin(iThetaY) 0 cos(iThetaY)) (0 sin(-iThetaX) cos(-iThetaX) )
      //
      //       (Xcos(iThetaY) + Zsin(iThetaY)cos(iThetaX) - Ysin(iThetaY)sin(iThetaX))
      // Out = (Ycos(iThetaX) + Zsin(iThetaX)                                        )
      //       (Zcos(iThetaX)cos(iThetaY) - Ysin(iThetaX)cos(iThetaY) - Xsin(iThetaY))

      const float lTemp = iZIn * iCosThetaX - iYIn * iSinThetaX;

      oXOut = lTemp * iSinThetaY + iXIn * iCosThetaY;
      oYOut = iZIn * iSinThetaX + iYIn * iCosThetaX;
      oZOut = lTemp * iCosThetaY - iXIn * iSinThetaY;
    }

    void Control::UndoRotateAboutXThenY(const float &iSinThetaX,
                                        const float &iCosThetaX,
                                        const float &iSinThetaY,
                                        const float &iCosThetaY,
                                        const float &iXIn,
                                        const float &iYIn,
                                        const float &iZIn,
                                        float &oXOut,
                                        float &oYOut,
                                        float &oZOut)
    {
      // Just do using the other Rotate function, swapping for x and y
      // (to swap order) and giving the inverse of the sine values
      // (to swap direction).
      RotateAboutXThenY(-iSinThetaY, iCosThetaY, -iSinThetaX, iCosThetaX, iYIn,
                        iXIn, iZIn, oYOut, oXOut, oZOut);
    }

    void Control::project(float p1[], float p2[])
    {
      float x1[3], x2[3];

      for (int l = 0; l < 3; l++)
      {
        x1[l] = p1[l] - mViewpoint.x[l];
      }

      UndoRotateAboutXThenY(mViewpoint.SinXRotation, mViewpoint.CosXRotation,
                            mViewpoint.SinYRotation, mViewpoint.CosYRotation,
                            x1[0], x1[1], x1[2], x2[0], x2[1], x2[2]);

      float temp = mViewpoint.dist / (p2[2] = -x2[2]);

      p2[0] = temp * x2[0];
      p2[1] = temp * x2[1];
    }

    void Control::SetProjection(const int &iPixels_x,
                                const int &iPixels_y,
                                const float &iLocal_ctr_x,
                                const float &iLocal_ctr_y,
                                const float &iLocal_ctr_z,
                                const float &iLongitude,
                                const float &iLatitude,
                                const float &iZoom)
    {
      mScreen.MaxXValue = (0.5 * vis->system_size) / iZoom;
      mScreen.MaxYValue = (0.5 * vis->system_size) / iZoom;

      mScreen.PixelsX = iPixels_x;
      mScreen.PixelsY = iPixels_y;

      // Convert to radians
      float temp = iLongitude * 0.01745329F;

      mViewpoint.SinYRotation = sinf(temp);
      mViewpoint.CosYRotation = cosf(temp);

      // Convert to radians
      temp = iLatitude * 0.01745329F;

      mViewpoint.SinXRotation = sinf(temp);
      mViewpoint.CosXRotation = cosf(temp);

      float rad = 5.F * vis->system_size;
      float dist = 0.5 * rad;

      temp = rad * mViewpoint.CosXRotation;

      mViewpoint.x[0] = temp * mViewpoint.SinYRotation + iLocal_ctr_x;
      mViewpoint.x[1] = rad * mViewpoint.SinXRotation + iLocal_ctr_y;
      mViewpoint.x[2] = temp * mViewpoint.CosYRotation + iLocal_ctr_z;

      mViewpoint.dist = dist;

      RotateAboutXThenY(mViewpoint.SinXRotation, mViewpoint.CosXRotation,
                        mViewpoint.SinYRotation, mViewpoint.CosYRotation,
                        mScreen.MaxXValue, 0.0F, 0.0F,
                        mScreen.UnitVectorProjectionX[0],
                        mScreen.UnitVectorProjectionX[1],
                        mScreen.UnitVectorProjectionX[2]);

      RotateAboutXThenY(mViewpoint.SinXRotation, mViewpoint.CosXRotation,
                        mViewpoint.SinYRotation, mViewpoint.CosYRotation, 0.0F,
                        mScreen.MaxYValue, 0.0F,
                        mScreen.UnitVectorProjectionY[0],
                        mScreen.UnitVectorProjectionY[1],
                        mScreen.UnitVectorProjectionY[2]);

      mScreen.ScaleX = (float) iPixels_x / (2.F * mScreen.MaxXValue);
      mScreen.ScaleY = (float) iPixels_y / (2.F * mScreen.MaxYValue);

      temp = dist / rad;

      mScreen.vtx[0] = (mViewpoint.x[0] + temp * (iLocal_ctr_x
          - mViewpoint.x[0])) - mScreen.UnitVectorProjectionX[0]
          - mScreen.UnitVectorProjectionY[0] - mViewpoint.x[0];
      mScreen.vtx[1] = (mViewpoint.x[1] + temp * (iLocal_ctr_y
          - mViewpoint.x[1])) - mScreen.UnitVectorProjectionX[1]
          - mScreen.UnitVectorProjectionY[1] - mViewpoint.x[1];
      mScreen.vtx[2] = (mViewpoint.x[2] + temp * (iLocal_ctr_z
          - mViewpoint.x[2])) - mScreen.UnitVectorProjectionX[2]
          - mScreen.UnitVectorProjectionY[2] - mViewpoint.x[2];

      mScreen.UnitVectorProjectionX[0] *= (2.F / (float) iPixels_x);
      mScreen.UnitVectorProjectionX[1] *= (2.F / (float) iPixels_x);
      mScreen.UnitVectorProjectionX[2] *= (2.F / (float) iPixels_x);

      mScreen.UnitVectorProjectionY[0] *= (2.F / (float) iPixels_y);
      mScreen.UnitVectorProjectionY[1] *= (2.F / (float) iPixels_y);
      mScreen.UnitVectorProjectionY[2] *= (2.F / (float) iPixels_y);
    }

    void Control::MergePixels(ColPixel *col_pixel1, ColPixel *col_pixel2)
    {
      // Merge raytracing data

      if (col_pixel1->i.isRt && col_pixel2->i.isRt)
      {
        // Both are ray-tracing
        col_pixel2->vel_r += col_pixel1->vel_r;
        col_pixel2->vel_g += col_pixel1->vel_g;
        col_pixel2->vel_b += col_pixel1->vel_b;

        if (mStressType != lb::ShearStress)
        {
          col_pixel2->stress_r += col_pixel1->stress_r;
          col_pixel2->stress_g += col_pixel1->stress_g;
          col_pixel2->stress_b += col_pixel1->stress_b;
        }

        col_pixel2->dt += col_pixel1->dt;

        if (col_pixel1->t < col_pixel2->t)
        {
          col_pixel2->t = col_pixel1->t;
          col_pixel2->density = col_pixel1->density;
          col_pixel2->stress = col_pixel1->stress;
        }

      }
      else if (col_pixel1->i.isRt && !col_pixel2->i.isRt)
      {
        // Only 1 is ray-tracing
        col_pixel2->vel_r = col_pixel1->vel_r;
        col_pixel2->vel_g = col_pixel1->vel_g;
        col_pixel2->vel_b = col_pixel1->vel_b;

        if (mStressType != lb::ShearStress)
        {
          col_pixel2->stress_r = col_pixel1->stress_r;
          col_pixel2->stress_g = col_pixel1->stress_g;
          col_pixel2->stress_b = col_pixel1->stress_b;
        }

        col_pixel2->t = col_pixel1->t;
        col_pixel2->dt = col_pixel1->dt;
        col_pixel2->density = col_pixel1->density;
        col_pixel2->stress = col_pixel1->stress;

        col_pixel2->i.isRt = true;
      }
      // Done merging ray-tracing

      // Now merge glyph data
      if (mStressType != lb::ShearStress && (mode == 0 || mode == 1))
      {

        if (col_pixel1->i.isGlyph)
        {
          col_pixel2->i.isGlyph = true;
        }

      }
      else
      {

#ifndef NO_STREAKLINES
        // merge streakline data
        if (col_pixel1->i.isStreakline && col_pixel2->i.isStreakline)
        {

          if (col_pixel1->particle_z < col_pixel2->particle_z)
          {
            col_pixel2->particle_z = col_pixel1->particle_z;
            col_pixel2->particle_vel = col_pixel1->particle_vel;
            col_pixel2->particle_inlet_id = col_pixel1->particle_inlet_id;
          }

        }
        else if (col_pixel1->i.isStreakline && !col_pixel2->i.isStreakline)
        {
          col_pixel2->particle_z = col_pixel1->particle_z;
          col_pixel2->particle_vel = col_pixel1->particle_vel;
          col_pixel2->particle_inlet_id = col_pixel1->particle_inlet_id;

          col_pixel2->i.isStreakline = true;
        }
#endif
      }

    }

    void Control::RegisterSite(int i,
                               float density,
                               float velocity,
                               float stress)
    {
      myRayTracer->UpdateClusterVoxel(i, density, velocity, stress);
    }

    void Control::writePixel(ColPixel *col_pixel_p)
    {
      int *col_pixel_id_p, i, j;

      i = col_pixel_p->i.i;
      j = col_pixel_p->i.j;

      col_pixel_id_p = &col_pixel_id[i * mScreen.PixelsY + j];

      if (*col_pixel_id_p != -1)
      {
        MergePixels(col_pixel_p, &col_pixel_send[*col_pixel_id_p]);

      }
      else
      { // col_pixel_id_p == -1

        if (col_pixels >= COLOURED_PIXELS_MAX)
        {
          return;
        }

        *col_pixel_id_p = col_pixels;

        memcpy(&col_pixel_send[col_pixels], col_pixel_p, sizeof(ColPixel));
        ++col_pixels;
      }

    }

    void Control::renderLine(float p1[], float p2[])
    {
      int pixels_x, pixels_y;
      int x, y;
      int x1, y1;
      int x2, y2;
      int dx, dy;
      int incE, incNE;
      int d;
      int m;

      ColPixel col_pixel;

      pixels_x = mScreen.PixelsX;
      pixels_y = mScreen.PixelsY;

      x1 = int(p1[0]);
      y1 = int(p1[1]);

      x2 = int(p2[0]);
      y2 = int(p2[1]);

      if (x2 < x1)
      {
        x = x1;
        y = y1;
        x1 = x2;
        y1 = y2;
        x2 = x;
        y2 = y;
      }

      x = x1;
      y = y1;

      if (y1 < y2)
      {
        dy = y2 - y1;
        m = 1;

      }
      else
      {
        dy = y1 - y2;
        m = -1;
      }

      dx = x2 - x1;

      if (dx > dy)
      {
        incE = dy;
        d = dy - dx;
        incNE = d;

        while (x <= x2)
        {
          if (! (x < 0 || x >= pixels_x || y < 0 || y >= pixels_y))
          {

            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            writePixel(&col_pixel);
          }

          if (d < 0)
          {
            d += incE;
          }
          else
          {
            d += incNE;
            y += m;
          }
          ++x;

        } // end while

      }
      else if (y1 < y2)
      {
        incE = dx;
        d = dx - dy;
        incNE = d;

        while (y <= y2)
        {
          if (! (x < 0 || x >= pixels_x || y < 0 || y >= pixels_y))
          {

            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            writePixel(&col_pixel);
          }

          if (d < 0)
          {
            d += incE;
          }
          else
          {
            d += incNE;
            ++x;
          }
          ++y;

        } // while

      }
      else
      {
        incE = dx;
        d = dx - dy;
        incNE = d;

        while (y >= y2)
        {
          if (! (x < 0 || x >= pixels_x || y < 0 || y >= pixels_y))
          {

            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            writePixel(&col_pixel);
          }

          if (d < 0)
          {
            d += incE;
          }
          else
          {
            d += incNE;
            ++x;
          }
          --y;

        } // while
      }
    }

    void Control::SetSomeParams(const float iBrightness,
                                const float iDensityThresholdMin,
                                const float iDensityThresholdMinMaxInv,
                                const float iVelocityThresholdMaxInv,
                                const float iStressThresholdMaxInv)
    {
      brightness = iBrightness;
      density_threshold_min = iDensityThresholdMin;

      density_threshold_minmax_inv = iDensityThresholdMinMaxInv;
      velocity_threshold_max_inv = iVelocityThresholdMaxInv;
      stress_threshold_max_inv = iStressThresholdMaxInv;
    }

    void Control::updateImageSize(int pixels_x, int pixels_y)
    {
      if (pixels_x * pixels_y > mScreen.PixelsX * mScreen.PixelsY)
      {
        pixels_max = pixels_x * pixels_y;
        col_pixel_id = (int *) realloc(col_pixel_id, sizeof(int) * pixels_max);
      }
      for (int i = 0; i < pixels_x * pixels_y; i++)
      {
        col_pixel_id[i] = -1;
      }
    }

    void Control::compositeImage(int recv_buffer_id,
                                 const topology::NetworkTopology * iNetTopology)
    {
      // here, the communications needed to composite the image are
      // handled through a binary tree pattern and parallel pairwise
      // blocking communications.

      MPI_Status status;

      int *col_pixel_id_p;
      int col_pixels_temp;
      int comm_inc, send_id, recv_id;
      int i, j;
      int m, n;

      ColPixel *col_pixel1, *col_pixel2;

#ifndef NEW_COMPOSITING
      memcpy(col_pixel_recv[recv_buffer_id], col_pixel_send, col_pixels
          * sizeof(ColPixel));
#else
      if (net->id != 0)
      {
        memcpy (col_pixel_recv[recv_buffer_id],
            col_pixel_send,
            col_pixels * sizeof(ColPixel));
      }
#endif

      comm_inc = 1;
      m = 1;

      while (m < iNetTopology->ProcessorCount)
      {
        m <<= 1;
#ifndef NEW_COMPOSITING
        int start_id = 0;
#else
        int start_id = 1;
#endif
        for (recv_id = start_id; recv_id < iNetTopology->ProcessorCount;)
        {
          send_id = recv_id + comm_inc;

          if (iNetTopology->LocalRank != recv_id && iNetTopology->LocalRank
              != send_id)
          {
            recv_id += comm_inc << 1;
            continue;
          }

          if (send_id >= iNetTopology->ProcessorCount || recv_id == send_id)
          {
            recv_id += comm_inc << 1;
            continue;
          }

          if (iNetTopology->LocalRank == send_id)
          {
            MPI_Send(&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);

            if (col_pixels > 0)
            {
              MPI_Send(col_pixel_send, col_pixels, ColPixel::getMpiType(),
                       recv_id, 20, MPI_COMM_WORLD);
            }

          }
          else
          {
            MPI_Recv(&col_pixels_temp, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
                     &status);

            if (col_pixels_temp > 0)
            {
              MPI_Recv(col_pixel_send, col_pixels_temp, ColPixel::getMpiType(),
                       send_id, 20, MPI_COMM_WORLD, &status);
            }

            for (n = 0; n < col_pixels_temp; n++)
            {
              col_pixel1 = &col_pixel_send[n];
              i = col_pixel1->i.i;
              j = col_pixel1->i.j;

              if (* (col_pixel_id_p = &col_pixel_id[i * mScreen.PixelsY + j])
                  == -1)
              {
                col_pixel2 = &col_pixel_recv[recv_buffer_id][*col_pixel_id_p
                    = col_pixels];

                memcpy(col_pixel2, col_pixel1, sizeof(ColPixel));
                ++col_pixels;

              }
              else
              {
                col_pixel2 = &col_pixel_recv[recv_buffer_id][*col_pixel_id_p];

                MergePixels(col_pixel1, col_pixel2);
              }

            }
          }
          if (m < iNetTopology->ProcessorCount && iNetTopology->LocalRank
              == recv_id)
          {
            memcpy(col_pixel_send, col_pixel_recv[recv_buffer_id], col_pixels
                * sizeof(ColPixel));
          }

          recv_id += comm_inc << 1;
        }
        comm_inc <<= 1;
      }
#ifdef NEW_COMPOSITING
      if (iNetTopology->LocalRank == 1)
      {
        MPI_Send (&col_pixels, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);

        if (col_pixels > 0)
        {
          MPI_Send (col_pixel_recv[recv_buffer_id], col_pixels, MPI_col_pixel_type,
              0, 20, MPI_COMM_WORLD);
        }

      }
      else if (iNetTopology->LocalRank == 0)
      {
        MPI_Recv (&col_pixels, 1, MPI_INT, 1, 20, MPI_COMM_WORLD, &status);

        if (col_pixels > 0)
        {
          MPI_Recv (col_pixel_recv[recv_buffer_id], col_pixels, MPI_col_pixel_type,
              1, 20, MPI_COMM_WORLD, &status);
        }

      }
#endif // NEW_COMPOSITING
    }

    void Control::render(int recv_buffer_id,
                         lb::GlobalLatticeData &iGlobLatDat,
                         const topology::NetworkTopology * iNetTopology)
    {
      if (mScreen.PixelsX * mScreen.PixelsY > pixels_max)
      {
        pixels_max = util::max(2 * pixels_max, mScreen.PixelsX
            * mScreen.PixelsY);

        col_pixel_id = (int *) realloc(col_pixel_id, sizeof(int) * pixels_max);
      }
      col_pixels = 0;

      myRayTracer->Render(mStressType);

      if (mode == 1)
      {
        myGlypher->render();
      }
#ifndef NO_STREAKLINES
      if (shouldDrawStreaklines
          && (mStressType == lb::ShearStress || mode == 2))
      {
        myStreaker->render(iGlobLatDat);
      }
#endif
      compositeImage(recv_buffer_id, iNetTopology);

      col_pixels_recv[recv_buffer_id] = col_pixels;
#ifndef NEW_COMPOSITING
      if (iNetTopology->IsCurrentProcTheIOProc())
      {
        return;
      }
#endif
      for (int m = 0; m < col_pixels_recv[recv_buffer_id]; m++)
      {
        col_pixel_id[col_pixel_send[m].i.i * mScreen.PixelsY
            + col_pixel_send[m].i.j] = -1;
      }
    }

    void Control::writeImage(int recv_buffer_id,
                             std::string image_file_name,
                             void(*ColourPalette)(float value, float col[]))
    {
      io::XdrFileWriter writer = io::XdrFileWriter(image_file_name);

      writer << mode;

      writer << physical_pressure_threshold_min
          << physical_pressure_threshold_max << physical_velocity_threshold_max
          << physical_stress_threshold_max;

      writer << mScreen.PixelsX << mScreen.PixelsY
          << col_pixels_recv[recv_buffer_id];

      for (int n = 0; n < col_pixels_recv[recv_buffer_id]; n++)
      {
        writer.writePixel(&col_pixel_recv[recv_buffer_id][n], ColourPalette,
                          mStressType);
      }
    }

    void Control::setMouseParams(double iPhysicalPressure,
                                 double iPhysicalStress)
    {
      mouse_pressure = iPhysicalPressure;
      mouse_stress = iPhysicalStress;
    }

    void Control::streaklines(int time_step,
                              int period,
                              lb::GlobalLatticeData &iGlobLatDat,
                              lb::LocalLatticeData &iLocalLatDat)
    {
      myStreaker ->StreakLines(time_step, period, iGlobLatDat, iLocalLatDat);
    }

    void Control::restart()
    {
      myStreaker->Restart();
    }

    Control::~Control()
    {
#ifndef NO_STREAKLINES
      delete myStreaker;
#endif

      delete vis;
      delete myGlypher;
      delete myRayTracer;

      delete[] col_pixel_recv[0];
      delete[] col_pixel_recv[1];

      delete[] col_pixel_id;
    }

  } // namespace vis
} // namespace hemelb
