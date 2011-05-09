#include "D3Q15.h"
#include <math.h>

namespace hemelb
{
  const int D3Q15::CX[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
  const int D3Q15::CY[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 };
  const int D3Q15::CZ[] = { 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1 };

  const int D3Q15::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13 };

  void D3Q15::CalculateDensityAndVelocity(const distribn_t f[],
                                          distribn_t &density,
                                          distribn_t &v_x,
                                          distribn_t &v_y,
                                          distribn_t &v_z)
  {
    v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
    v_y = f[3] + (f[12] + f[14]);
    v_z = f[5] + f[10];

    density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + v_x + v_y + v_z;

    v_x -= (f[2] + f[8] + f[10] + (f[12] + f[14]));
    v_y += (f[7] + f[9]) - ( (f[4] + f[8] + f[10] + (f[11] + f[13])));
    v_z += f[7] + f[11] + f[14] - ( ( (f[6] + f[8]) + f[9] + f[12] + f[13]));
  }

  void D3Q15::CalculateFeq(const distribn_t &density,
                           const distribn_t &v_x,
                           const distribn_t &v_y,
                           const distribn_t &v_z,
                           distribn_t f_eq[])
  {
    distribn_t density_1 = 1. / density;
    distribn_t v_xx = v_x * v_x;
    distribn_t v_yy = v_y * v_y;
    distribn_t v_zz = v_z * v_z;

    f_eq[0] = (2.0 / 9.0) * density - (1.0 / 3.0) * ( (v_xx + v_yy + v_zz) * density_1);

    distribn_t temp1 = (1.0 / 9.0) * density - (1.0 / 6.0) * ( (v_xx + v_yy + v_zz) * density_1);

    f_eq[1] = (temp1 + (0.5 * density_1) * v_xx) + ( (1.0 / 3.0) * v_x); // (+1, 0, 0)
    f_eq[2] = (temp1 + (0.5 * density_1) * v_xx) - ( (1.0 / 3.0) * v_x); // (+1, 0, 0)

    f_eq[3] = (temp1 + (0.5 * density_1) * v_yy) + ( (1.0 / 3.0) * v_y); // (0, +1, 0)
    f_eq[4] = (temp1 + (0.5 * density_1) * v_yy) - ( (1.0 / 3.0) * v_y); // (0, -1, 0)

    f_eq[5] = (temp1 + (0.5 * density_1) * v_zz) + ( (1.0 / 3.0) * v_z); // (0, 0, +1)
    f_eq[6] = (temp1 + (0.5 * density_1) * v_zz) - ( (1.0 / 3.0) * v_z); // (0, 0, -1)

    temp1 *= (1.0 / 8.0);

    distribn_t temp2 = (v_x + v_y) + v_z;

    f_eq[7] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, +1, +1)
    f_eq[8] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, -1, -1)

    temp2 = (v_x + v_y) - v_z;

    f_eq[9] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, +1, -1)
    f_eq[10] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, -1, +1)

    temp2 = (v_x - v_y) + v_z;

    f_eq[11] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, -1, +1)
    f_eq[12] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, +1, -1)

    temp2 = (v_x - v_y) - v_z;

    f_eq[13] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, -1, -1)
    f_eq[14] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, +1, +1)
  }

  // Calculate density, velocity and the equilibrium distribution
  // functions according to the D3Q15 model.  The calculated v_x, v_y
  // and v_z are actually density * velocity, because we are using the
  // compressible model.
  // TODO Alter this to make it return actual velocity.
  void D3Q15::CalculateDensityVelocityFEq(const distribn_t f[],
                                          distribn_t &density,
                                          distribn_t &v_x,
                                          distribn_t &v_y,
                                          distribn_t &v_z,
                                          distribn_t f_eq[])
  {
    v_x = f[1] + (f[7] + f[9]) + (f[11] + f[13]);
    v_y = f[3] + (f[12] + f[14]);
    v_z = f[5] + f[10];

    density = f[0] + (f[2] + f[4]) + (f[6] + f[8]) + v_x + v_y + v_z;

    v_x -= f[2] + (f[8] + f[10]) + (f[12] + f[14]);
    v_y += (f[7] + f[9]) - (f[4] + (f[8] + f[10]) + (f[11] + f[13]));
    v_z += f[7] + f[11] + f[14] - ( (f[6] + f[8]) + f[9] + f[12] + f[13]);

    distribn_t density_1 = 1. / density;

    distribn_t v_xx = v_x * v_x;
    distribn_t v_yy = v_y * v_y;
    distribn_t v_zz = v_z * v_z;

    f_eq[0] = (2.0 / 9.0) * density - (1.0 / 3.0) * ( (v_xx + v_yy + v_zz) * density_1);

    distribn_t temp1 = (1.0 / 9.0) * density - (1.0 / 6.0) * ( (v_xx + v_yy + v_zz) * density_1);

    f_eq[1] = (temp1 + (0.5 * density_1) * v_xx) + ( (1.0 / 3.0) * v_x); // (+1, 0, 0)
    f_eq[2] = (temp1 + (0.5 * density_1) * v_xx) - ( (1.0 / 3.0) * v_x); // (+1, 0, 0)

    f_eq[3] = (temp1 + (0.5 * density_1) * v_yy) + ( (1.0 / 3.0) * v_y); // (0, +1, 0)
    f_eq[4] = (temp1 + (0.5 * density_1) * v_yy) - ( (1.0 / 3.0) * v_y); // (0, -1, 0)

    f_eq[5] = (temp1 + (0.5 * density_1) * v_zz) + ( (1.0 / 3.0) * v_z); // (0, 0, +1)
    f_eq[6] = (temp1 + (0.5 * density_1) * v_zz) - ( (1.0 / 3.0) * v_z); // (0, 0, -1)

    temp1 *= (1.0 / 8.0);

    double temp2 = (v_x + v_y) + v_z;

    f_eq[7] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, +1, +1)
    f_eq[8] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, -1, -1)

    temp2 = (v_x + v_y) - v_z;

    f_eq[9] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, +1, -1)
    f_eq[10] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, -1, +1)

    temp2 = (v_x - v_y) + v_z;

    f_eq[11] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, -1, +1)
    f_eq[12] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, +1, -1)

    temp2 = (v_x - v_y) - v_z;

    f_eq[13] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) + ( (1.0 / 24.0) * temp2); // (+1, -1, -1)
    f_eq[14] = (temp1 + ( (1.0 / 16.0) * density_1) * temp2 * temp2) - ( (1.0 / 24.0) * temp2); // (-1, +1, +1)
  }

  // von Mises stress computation given the non-equilibrium distribution functions.
  void D3Q15::CalculateVonMisesStress(const distribn_t f[],
                                      distribn_t &stress,
                                      const double iStressParameter)
  {
    distribn_t sigma_xx_yy = (f[1] + f[2]) - (f[3] + f[4]);
    distribn_t sigma_yy_zz = (f[3] + f[4]) - (f[5] + f[6]);
    distribn_t sigma_xx_zz = (f[1] + f[2]) - (f[5] + f[6]);

    distribn_t sigma_xy = (f[7] + f[8]) + (f[9] + f[10]) - (f[11] + f[12]) - (f[13] + f[14]);
    distribn_t sigma_xz = (f[7] + f[8]) - (f[9] + f[10]) + (f[11] + f[12]) - (f[13] + f[14]);
    distribn_t sigma_yz = (f[7] + f[8]) - (f[9] + f[10]) - (f[11] + f[12]) + (f[13] + f[14]);

    distribn_t a = sigma_xx_yy * sigma_xx_yy + sigma_yy_zz * sigma_yy_zz + sigma_xx_zz
        * sigma_xx_zz;
    distribn_t b = sigma_xy * sigma_xy + sigma_xz * sigma_xz + sigma_yz * sigma_yz;

    stress = iStressParameter * sqrt(a + 6.0 * b);
  }

  // The magnitude of the tangential component of the shear stress acting on the
  // wall.
  void D3Q15::CalculateShearStress(const distribn_t &density,
                                   const distribn_t f[],
                                   const double nor[],
                                   distribn_t &stress,
                                   const double &iStressParameter)
  {
    distribn_t sigma[9]; // stress tensor;
    // sigma_ij is the force
    // per unit area in
    // direction i on the
    // plane with the normal
    // in direction j
    distribn_t stress_vector[] = { 0.0, 0.0, 0.0 }; // Force per unit area in
    // direction i on the
    // plane perpendicular to
    // the surface normal
    distribn_t square_stress_vector = 0.0;
    distribn_t normal_stress = 0.0; // Magnitude of force per
    // unit area normal to the
    // surface
    unsigned int i, j, l;

    distribn_t temp = iStressParameter * (-sqrt(2.0));

    const int *Cs[3] = { CX, CY, CZ };

    for (i = 0; i < 3; i++)
    {
      for (j = 0; j <= i; j++)
      {
        sigma[i * 3 + j] = 0.0;

        for (l = 0; l < NUMVECTORS; l++)
        {
          sigma[i * 3 + j] += f[l] * (Cs[i][l] * Cs[j][l]);
        }
        sigma[i * 3 + j] *= temp;
      }
    }
    for (i = 0; i < 3; i++)
      for (j = 0; j < i; j++)
        sigma[j * 3 + i] = sigma[i * 3 + j];

    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
        stress_vector[i] += sigma[i * 3 + j] * nor[j];

      square_stress_vector += stress_vector[i] * stress_vector[i];
      normal_stress += stress_vector[i] * nor[i];
    }
    // shear_stress^2 + normal_stress^2 = stress_vector^2
    stress = sqrt(square_stress_vector - normal_stress * normal_stress);
  }
}
