#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits>

#include "vis/rayTracer/HSLToRGBConverter.h"
//#include "HSLToRGBConverter.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      //See http://en.wikipedia.org/wiki/HSL_and_HSV for the
      //formulae used
      //TODO: Possibly replace truncation with rounding?
      void HSLToRGBConverter::Convert(float iHue,
                                      float iSaturation,
                                      float iLightness,
                                      unsigned char oRGBColour[3])
      {
        //All values are stored as 32-bit unsigned integers
        //stored within 16-bits, to allow multiplication between
        //two values 

        //All values are scaled up to be between 0 and the
        //maximum value of a 16-bit interger

        //Scales a value between 0 and 360 to a range between 0
        //and the maximum value of a 16-bit interger
        static const uint32_t DegreesScalar = std::numeric_limits<uint16_t>::max() / 360;

        //Scales a value between 0 and 1 to a range between 0
        //and the maximum value of a 16-bit interger
        static const uint32_t OtherScalar = std::numeric_limits<uint16_t>::max();

        //Cast the inputs accodingly
        uint32_t lHue = (uint32_t)(iHue * (float) (DegreesScalar));
        uint32_t lSaturation = (uint32_t)(iSaturation * (float) (OtherScalar));
        uint32_t lLightness = (uint32_t)(iLightness * (float) (OtherScalar));

        //Calculate the Chroma - a division by OtherScalar is 
        //required for the Chroma to remain between 0
        //and the maximum value of a 16-bit interger 
        int32_t lTemp = 2 * (int32_t) (lLightness) - (int32_t) (OtherScalar);
        uint32_t lChroma = ( (OtherScalar - abs(lTemp)) * lSaturation) / OtherScalar;

        uint32_t lHuePrime = lHue / 60;

        lTemp = (int) (lHuePrime % (2 * DegreesScalar)) - (int) (DegreesScalar);

        //Calculate the "Intermediate" - a division by 
        //DegreesScalar is required for the Chroma to 
        //remain between 0 and the maximum value of a 16-bit interger 
        uint32_t lIntermediate = (lChroma * (DegreesScalar - abs(lTemp))) / DegreesScalar;

        uint32_t lRed;
        uint32_t lGreen;
        uint32_t lBlue;

        //Map the hue to value to six cases 
        uint32_t lHueInt = lHuePrime / DegreesScalar;
        switch (lHueInt)
        {
          case 0:
            lRed = lChroma;
            lGreen = lIntermediate;
            lBlue = 0.0F;
            break;

          case 1:
            lRed = lIntermediate;
            lGreen = lChroma;
            lBlue = 0.0F;
            break;

          case 2:
            lRed = 0.0F;
            lGreen = lChroma;
            lBlue = lIntermediate;
            break;

          case 3:
            lRed = 0.0F;
            lGreen = lIntermediate;
            lBlue = lChroma;
            break;

          case 4:
            lRed = lIntermediate;
            lGreen = 0.0F;
            lBlue = lChroma;
            break;

            // lHueInt is guaranteed to be 0-5, but having a default case rather than 'case 5:'
            // prevents a compiler warning (that lRed, lGreen and lBlue may not have been
            // assigned)
          default:
            lRed = lChroma;
            lGreen = 0.0F;
            lBlue = lIntermediate;
            break;
        }

        // A value to divide the results by to map them 
        // between 0 and the maximum value of an unsigned
        //char
        static const uint32_t CharScalar = OtherScalar / std::numeric_limits<unsigned char>::max();

        int32_t lMatcher = (int) (lLightness) - (int) (lChroma) / 2;

        oRGBColour[0] = (unsigned char) ( (lRed + lMatcher) / CharScalar);
        oRGBColour[1] = (unsigned char) ( (lGreen + lMatcher) / CharScalar);
        oRGBColour[2] = (unsigned char) ( (lBlue + lMatcher) / CharScalar);
      }
    }
  }
}
