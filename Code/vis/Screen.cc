#include <stdlib.h>

#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {

    // TODO This is probably going to have to be cleverly redesigned. We need to pass the images around over several iterations without
    // interference between the steering and written-to-disk images.

    Screen::Screen()
    {
    }

    Screen::~Screen()
    {
    }

    /**
     * Add a pixel to the screen.
     *
     * @param newPixel The new pixel to be added
     * @param iStressType The stress type of the visualisation
     * @param mode Controls what aspects of the visualisation to display.
     */
    void Screen::AddPixel(const ColPixel* newPixel, const VisSettings* visSettings)
    {
      localPixels.AddPixel(newPixel, visSettings);
    }

    /**
     * Render a line between two points on the screen.
     *
     * @param endPoint1
     * @param endPoint2
     * @param iStressType
     * @param mode
     */
    void Screen::RenderLine(const float endPoint1[3],
                            const float endPoint2[3],
                            const VisSettings* visSettings)
    {
      // Store end points of the line and 'current' point (x and y).
      int x = int (endPoint1[0]);
      int y = int (endPoint1[1]);

      int x1, y1, x2, y2;
      if (endPoint2[0] < endPoint1[0])
      {
        x1 = int (endPoint2[0]);
        y1 = int (endPoint2[1]);
        x2 = x;
        y2 = y;
      }
      else
      {
        x1 = x;
        y1 = y;
        x2 = int (endPoint2[0]);
        y2 = int (endPoint2[1]);
      }

      // Set dx with the difference between x-values at the endpoints.
      int dx = x2 - x1;

      // Set up the iteration in general terms.

      if (y2 <= y1)
      {
        int temp = y2;
        y2 = y;
        y = temp;
      }

      int dy = y2 - y;

      if (dx > dy)
      {
        RenderLineHelper<true> (x, y, dy, dy - dx, x2, visSettings);
      }
      else
      {
        RenderLineHelper<false> (x, y, dx, dx - dy, y2, visSettings);
      }
    }

    /**
     * Helper function for rendering a line. The template parameter indicates whether the end of
     * the line will be limited by the x or y dimension.
     *
     * @param x
     * @param y
     * @param incE
     * @param incNE
     * @param limit
     * @param stressType
     * @param mode
     */
    template<bool xLimited>
    void Screen::RenderLineHelper(int x,
                                  int y,
                                  int incE,
                                  int incNE,
                                  int limit,
                                  const VisSettings* visSettings)
    {
      int d = incNE;

      while ( (xLimited && x <= limit) || (!xLimited && y <= limit))
      {
        if (x >= 0 && x < (int) localPixels.GetPixelsX() && y >= 0 && y
            < (int) localPixels.GetPixelsY())
        {
          ColPixel col_pixel(x, y);
          AddPixel(&col_pixel, visSettings);
        }

        if (d < 0)
        {
          d += incE;
          if (xLimited)
          {
            ++x;
          }
          else
          {
            ++y;
          }
        }
        else
        {
          d += incNE;
          ++y;
          ++x;
        }
      } // end while
    }

    void Screen::Set(float maxX,
                     float maxY,
                     int pixelsX,
                     int pixelsY,
                     float rad,
                     const Viewpoint* viewpoint)
    {
      MaxXValue = maxX;
      MaxYValue = maxX;

      viewpoint->RotateToViewpoint(MaxXValue, 0.0F, 0.0F, UnitVectorProjectionX);
      viewpoint->RotateToViewpoint(0.0F, MaxYValue, 0.0F, UnitVectorProjectionY);

      if (pixelsX * pixelsY <= (int) COLOURED_PIXELS_MAX)
      {
        localPixels.SetSize(pixelsX, pixelsY);
      }

      ScaleX = (float) localPixels.GetPixelsX() / (2.F * MaxXValue);
      ScaleY = (float) localPixels.GetPixelsY() / (2.F * MaxYValue);

      float radVector[3];
      viewpoint->RotateToViewpoint(0.F, 0.F, -rad, radVector);

      for (int ii = 0; ii < 3; ++ii)
      {
        vtx[ii] = (0.5F * radVector[ii]) - UnitVectorProjectionX[ii] - UnitVectorProjectionY[ii];
      }

      UnitVectorProjectionX[0] *= (2.F / (float) localPixels.GetPixelsX());
      UnitVectorProjectionX[1] *= (2.F / (float) localPixels.GetPixelsX());
      UnitVectorProjectionX[2] *= (2.F / (float) localPixels.GetPixelsX());

      UnitVectorProjectionY[0] *= (2.F / (float) localPixels.GetPixelsY());
      UnitVectorProjectionY[1] *= (2.F / (float) localPixels.GetPixelsY());
      UnitVectorProjectionY[2] *= (2.F / (float) localPixels.GetPixelsY());
    }

    void Screen::Resize(unsigned int newPixelsX, unsigned int newPixelsY)
    {
      localPixels.SetSize(newPixelsX, newPixelsY);
    }

    void Screen::Reset()
    {
      localPixels.pixelCount = 0;
    }

    void Screen::CompositeImage(const VisSettings* visSettings)
    {
      // Status object for MPI comms.
      MPI_Status status;

      /*
       * We do several iterations.
       *
       * On the first, every even proc passes data to the odd proc below, where it is merged.
       * On the second, the difference is two, so proc 3 passes to 1, 7 to 5, 11 to 9 etc.
       * On the third the differenec is four, so proc 5 passes to 1, 13 to 9 etc.
       * .
       * .
       * .
       *
       * This continues until all data is passed back to processor one, which passes it to proc 0.
       */

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      // Start with a difference in rank of 1, doubling every time.
      for (proc_t deltaRank = 1; deltaRank < netTop->GetProcessorCount(); deltaRank <<= 1)
      {
        // The receiving proc is all the ranks that are 1 modulo (deltaRank * 2)
        for (proc_t receivingProc = 1; receivingProc < (netTop->GetProcessorCount() - deltaRank); receivingProc
            += deltaRank << 1)
        {
          proc_t sendingProc = receivingProc + deltaRank;

          // If we're the sending proc, do the send.
          if (netTop->GetLocalRank() == sendingProc)
          {
            MPI_Send(&localPixels.pixelCount,
                     1,
                     MpiDataType(localPixels.pixelCount),
                     receivingProc,
                     20,
                     MPI_COMM_WORLD);

            if (localPixels.pixelCount > 0)
            {
              MPI_Send(localPixels.pixels,
                       localPixels.pixelCount,
                       MpiDataType<ColPixel> (),
                       receivingProc,
                       20,
                       MPI_COMM_WORLD);
            }
          }

          // If we're the receiving proc, receive.
          else if (netTop->GetLocalRank() == receivingProc)
          {
            MPI_Recv(&compositingBuffer.pixelCount,
                     1,
                     MpiDataType(compositingBuffer.pixelCount),
                     sendingProc,
                     20,
                     MPI_COMM_WORLD,
                     &status);

            if (compositingBuffer.pixelCount > 0)
            {
              MPI_Recv(compositingBuffer.pixels,
                       compositingBuffer.pixelCount,
                       MpiDataType<ColPixel> (),
                       sendingProc,
                       20,
                       MPI_COMM_WORLD,
                       &status);
            }

            localPixels.AddPixels(compositingBuffer.pixels,
                                  compositingBuffer.pixelCount,
                                  visSettings);
          }
        }
      }

      // Send the final image from proc 1 to 0.
      if (netTop->GetLocalRank() == 1)
      {
        MPI_Send(&localPixels.pixelCount,
                 1,
                 MpiDataType(localPixels.pixelCount),
                 0,
                 20,
                 MPI_COMM_WORLD);

        if (localPixels.pixelCount > 0)
        {
          MPI_Send(localPixels.pixels,
                   localPixels.pixelCount,
                   MpiDataType<ColPixel> (),
                   0,
                   20,
                   MPI_COMM_WORLD);
        }

      }
      // Receive the final image on proc 0.
      else if (netTop->GetLocalRank() == 0)
      {
        MPI_Recv(&localPixels.pixelCount,
                 1,
                 MpiDataType(localPixels.pixelCount),
                 1,
                 20,
                 MPI_COMM_WORLD,
                 &status);

        if (localPixels.pixelCount > 0)
        {
          MPI_Recv(localPixels.pixels,
                   localPixels.pixelCount,
                   MpiDataType<ColPixel> (),
                   1,
                   20,
                   MPI_COMM_WORLD,
                   &status);
        }
      }

      pixelCountInBuffer = localPixels.pixelCount;

      for (unsigned int m = 0; m < pixelCountInBuffer; m++)
      {
        localPixels.pixelId[localPixels.pixels[m].GetI() * GetPixelsY()
            + localPixels.pixels[m].GetJ()] = -1;
      }
    }

    bool Screen::MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress)
    {
      for (unsigned int i = 0; i < pixelCountInBuffer; i++)
      {
        if (localPixels.pixels[i].IsRT() && int (localPixels.pixels[i].GetI()) == mouseX
            && int (localPixels.pixels[i].GetJ()) == mouseY)
        {
          *density = localPixels.pixels[i].GetDensity();
          *stress = localPixels.pixels[i].GetStress();

          return true;
        }
      }

      return false;
    }

    void Screen::WritePixelCount(io::Writer* writer)
    {
      writer->operator <<(GetPixelsX());
      writer->operator <<(GetPixelsY());
      writer->operator <<(pixelCountInBuffer);
    }

    void Screen::WritePixels(const DomainStats* domainStats,
                             const VisSettings* visSettings,
                             io::Writer* writer)
    {
      int index;
      int pix_data[3];
      unsigned char rgb_data[12];
      int bits_per_char = sizeof(char) * 8;

      for (unsigned int i = 0; i < pixelCountInBuffer; i++)
      {
        ColPixel& col_pixel = localPixels.pixels[i];
        // Use a ray-tracer function to get the necessary pixel data.
        col_pixel.rawWritePixel(&index, rgb_data, domainStats, visSettings);

        *writer << index;

        pix_data[0] = (rgb_data[0] << (3 * bits_per_char)) + (rgb_data[1] << (2 * bits_per_char))
            + (rgb_data[2] << bits_per_char) + rgb_data[3];

        pix_data[1] = (rgb_data[4] << (3 * bits_per_char)) + (rgb_data[5] << (2 * bits_per_char))
            + (rgb_data[6] << bits_per_char) + rgb_data[7];

        pix_data[2] = (rgb_data[8] << (3 * bits_per_char)) + (rgb_data[9] << (2 * bits_per_char))
            + (rgb_data[10] << bits_per_char) + rgb_data[11];

        for (int i = 0; i < 3; i++)
        {
          *writer << pix_data[i];
        }
        *writer << io::Writer::eol;
      }
    }

    unsigned int Screen::GetPixelCount() const
    {
      return pixelCountInBuffer;
    }

    const float* Screen::GetVtx() const
    {
      return vtx;
    }
    const float* Screen::GetUnitVectorProjectionX() const
    {
      return UnitVectorProjectionX;
    }
    const float* Screen::GetUnitVectorProjectionY() const
    {
      return UnitVectorProjectionY;
    }
    int Screen::GetPixelsX() const
    {
      return localPixels.GetPixelsX();
    }
    int Screen::GetPixelsY() const
    {
      return localPixels.GetPixelsY();
    }
  }
}
