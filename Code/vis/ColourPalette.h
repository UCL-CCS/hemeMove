#ifndef HEME_VIS_COLOURPALETTE_H
#define HEME_VIS_COLOURPALETTE_H

namespace heme
{
  namespace vis
  {
  
    // Type for conversion functions.
    typedef void (ColourPaletteFunction)(float, float *);

    class ColourPalette {
    public:
      // Function to populate a RGB colour (col) with values dependant
      //on the value of t.
      static ColourPaletteFunction pickColour;
    };
  
  }
}

#endif // HEME_VIS_COLOURPALETTE_H
