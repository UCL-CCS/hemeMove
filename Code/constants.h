#ifndef HEMELB_CONSTANTS_H
#define HEMELB_CONSTANTS_H

#include <limits>
#include <stdint.h>
#include "units.h"
//#include "mpiInclude.h"

namespace hemelb
{

  const unsigned int COLLISION_TYPES = 6;
  const double PI = 3.14159265358979323846264338327950288;
  const double DEG_TO_RAD = (PI / 180.0);
  // TODO this was used for a convergence test - we could reinstate that at some point.
  const double EPSILON = 1.0e-30;

  const double REFERENCE_PRESSURE_mmHg = 80.0;
  const double mmHg_TO_PASCAL = 133.3223874;
  const double BLOOD_DENSITY_Kg_per_m3 = 1000.0;
  const double BLOOD_VISCOSITY_Pa_s = 0.004;
  const double PULSATILE_PERIOD_s = 60.0 / 70.0;

  // These constants are used to pack a lot of stuff into a 32 bit int.
  // They are used in the setup tool and must be consistent.

  /* This is the number of boundary types. It was 4, but the
   * "CHARACTERISTIC_BOUNDARY" type is never used and I don't know what it is
   * meant to be. It is also not used in the setup tool, so we will drop it,
   * setting BOUNDARIES to 3
   */
  const sitedata_t BOUNDARIES = 3U;
  const sitedata_t INLET_BOUNDARY = 0U;
  const sitedata_t OUTLET_BOUNDARY = 1U;
  const sitedata_t WALL_BOUNDARY = 2U;
  // const unsigned int CHARACTERISTIC_BOUNDARY = 3U;

  const sitedata_t SITE_TYPE_BITS = 2U;
  const sitedata_t BOUNDARY_CONFIG_BITS = 26U;
  const sitedata_t BOUNDARY_DIR_BITS = 4U;
  const sitedata_t BOUNDARY_ID_BITS = 10U;

  const sitedata_t BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS;
  const sitedata_t BOUNDARY_DIR_SHIFT = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
  const sitedata_t BOUNDARY_ID_SHIFT = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

  // Comments show the bit patterns.
  const sitedata_t SITE_TYPE_MASK = ( (1 << SITE_TYPE_BITS) - 1);
  // 0000 0000  0000 0000  0000 0000  0000 0011
  // These give the *_TYPE in geometry/LatticeData.h

  const sitedata_t BOUNDARY_CONFIG_MASK = ( (1 << BOUNDARY_CONFIG_BITS) - 1)
      << BOUNDARY_CONFIG_SHIFT;
  // 0000 0000  0000 0000  1111 1111  1111 1100
  // These bits are set if the lattice vector they correspond to takes one to a solid site
  // The following hex digits give the index into LatticeSite.neighbours
  // ---- ----  ---- ----  DCBA 9876  5432 10--

  const sitedata_t BOUNDARY_DIR_MASK = ( (1 << BOUNDARY_DIR_BITS) - 1) << BOUNDARY_DIR_SHIFT;
  // 0000 0000  0000 1111  0000 0000  0000 0000
  // No idea what these represent. As far as I can tell, they're unused.

  const sitedata_t BOUNDARY_ID_MASK = ( ((sitedata_t) 1 << BOUNDARY_ID_BITS) - 1)
      << BOUNDARY_ID_SHIFT;
  // 0011 1111  1111 0000  0000 0000  0000 0000
  // These bits together give the index of the inlet/outlet/wall in the output XML file

  const sitedata_t PRESSURE_EDGE_MASK = (sitedata_t) 1
      << (BOUNDARY_ID_BITS + BOUNDARY_ID_SHIFT + 1);
  // 1000 0000  0000 0000  0000 0000  0000 0000

  const unsigned int FLUID = 1U;
  const unsigned int INLET = 2U;
  const unsigned int OUTLET = 4U;
  const unsigned int EDGE = 8U;

  // square of the speed of sound
  const double Cs2 = 1.0 / 3.0;

  // TODO almost certainly filth.
  const distribn_t NO_VALUE = std::numeric_limits<distribn_t>::max();
  const int BIG_NUMBER2 = 1 << 30;
  const unsigned BIG_NUMBER3 = 1 << 31;
}

#endif //HEMELB_CONSTANTS_H
