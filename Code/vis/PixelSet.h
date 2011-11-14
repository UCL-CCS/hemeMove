#ifndef HEMELB_VIS_PIXELSET_H
#define HEMELB_VIS_PIXELSET_H

#include <vector>
#include <map>

#include "log/Logger.h"
#include "mpiInclude.h"
#include "net/net.h"
#include "vis/BasicPixel.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    /**
     * Base pixel set implementation, including the functionality that allows storing of pixels
     * in a vector for speedy MPI usage, but also storing a map of pixels -> Pixel* for efficient
     * lookup in O(N) time.
     */
    template<typename PixelType>
    class PixelSet
    {
      public:

        PixelSet()
        {
          inUse = false;
          count = 0;
        }

        ~PixelSet()
        {
        }

        void Combine(const PixelSet<PixelType> &other)
        {
          typename std::vector<PixelType>::const_iterator otherPixels = other.pixels.begin();

          while (otherPixels != other.pixels.end())
          {
            AddPixel(*otherPixels);
            otherPixels++;
          }
        }

        void AddPixel(const PixelType& newPixel)
        {
          BasicPixel location = BasicPixel(newPixel.GetI(), newPixel.GetJ());

          if (pixelLookup.count(location) > 0)
          {
            pixels[pixelLookup[location]].Combine(newPixel);
          }
          else
          {
            pixels.push_back(PixelType(newPixel));
            pixelLookup.insert(std::pair<BasicPixel, unsigned int>(location,
                                                                   (unsigned int) (pixels.size()
                                                                       - 1)));
          }
        }

        bool IsInUse() const
        {
          return inUse;
        }

        void SetInUse()
        {
          inUse = true;
        }

        void Release()
        {
          inUse = false;
        }

        void SendQuantity(net::Net* net, proc_t destination)
        {
          count = (int) pixels.size();
          log::Logger::Log<log::Debug, log::OnePerCore>("Sending pixel count of %i", count);
          net->RequestSend(&count, 1, destination);
        }

        void ReceiveQuantity(net::Net* net, proc_t source)
        {
          net->RequestReceive(&count, 1, source);
        }

        void SendPixels(net::Net* net, proc_t destination)
        {
          if (pixels.size() > 0)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Sending %i pixels to proc %i",
                                                          (int) pixels.size(),
                                                          (int) destination);
            net->RequestSend(&pixels[0], (int) pixels.size(), destination);
          }
        }

        void ReceivePixels(net::Net* net, proc_t source)
        {
          if (count > 0)
          {
            // First make sure the vector will be large enough to hold the incoming pixels.
            log::Logger::Log<log::Debug, log::OnePerCore>("Receiving %i pixels from proc %i",
                                                          count,
                                                          (int) source);
            pixels.resize(count);
            net->RequestReceive(&pixels[0], count, source);
          }
        }

        size_t GetPixelCount() const
        {
          return pixels.size();
        }

        const std::vector<PixelType>& GetPixels() const
        {
          return pixels;
        }

        void Clear()
        {
          pixelLookup.clear();
          pixels.clear();
        }

      private:
        std::map<BasicPixel, unsigned int> pixelLookup;
        std::vector<PixelType> pixels;
        int count;
        bool inUse;
    };
  }
}

#endif // HEMELB_VIS_PIXELSET_H
