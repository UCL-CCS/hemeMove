// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetFileVelocity.h"

#include <algorithm>
#include <fstream>
#include "hassert.h"
#include "log/Logger.h"
#include "util/numerical.h"
#include "configuration/SimConfig.h"
#include <cmath>

namespace hemelb::lb
{
    InOutLetFileVelocity::InOutLetFileVelocity() :
            units(nullptr)
    {
    }

    InOutLet* InOutLetFileVelocity::clone() const
    {
        InOutLet* copy = new InOutLetFileVelocity(*this);
        return copy;
    }

    template<typename T>
    T LinearInterpolate(std::vector<T> const& xVector, std::vector<T> const& yVector, T targetX)
    {
        int lowerIndex = 0;

        if (log::Logger::ShouldDisplay<log::Debug>())
        {
            if (targetX < xVector[0] || targetX > xVector[xVector.size() - 1])
            {
                log::Logger::Log<log::Warning, log::OnePerCore>("Linear Interpolation beyond bounds: %f is not between %f and %f",
                                                                targetX,
                                                                xVector[0],
                                                                xVector[xVector.size() - 1]);
            }
        }

        while (! (targetX >= xVector[lowerIndex] && targetX <= xVector[lowerIndex + 1]))
        {
            lowerIndex++;
        }

        // If discontinuities are present in the trace the correct behaviour is ill-defined.
        if (log::Logger::ShouldDisplay<log::Debug>())
        {
            if (xVector[lowerIndex] == xVector[lowerIndex + 1])
            {
                log::Logger::Log<log::Warning, log::OnePerCore>("Multiple points for same x value in LinearInterpolate: ");
                log::Logger::Log<log::Warning, log::OnePerCore>("(%f, %f) and (%f, %f). Division by zero!",
                                                                xVector[lowerIndex],
                                                                yVector[lowerIndex],
                                                                xVector[lowerIndex + 1],
                                                                yVector[lowerIndex + 1]);
            }
        }

        // Linear interpolation of function f(x) between two points A and B
        // f(A) + (fraction along x axis between A and B) * (f(B) - f(A))
        return std::lerp(yVector[lowerIndex], yVector[lowerIndex + 1],
                         (targetX - xVector[lowerIndex]) / (xVector[lowerIndex + 1] - xVector[lowerIndex]));
    }

    void InOutLetFileVelocity::CalculateTable(LatticeTimeStep startTS, LatticeTimeStep endTS, PhysicalTime timeStepLength)
    {
        // First read in values from file
        // Used to be complex code here to keep a vector unique, but this is just achieved by using a map.
        std::map<PhysicalTime, PhysicalSpeed> timeValuePairs;

        double timeTemp, valueTemp;

        if (!std::filesystem::exists(velocityFilePath))
            throw Exception() << "File does not exist: " << velocityFilePath;

        std::ifstream datafile(velocityFilePath);
        log::Logger::Log<log::Debug, log::OnePerCore>("Reading iolet values from file:");
        while (datafile.good())
        {
          datafile >> timeTemp >> valueTemp;
          log::Logger::Log<log::Trace, log::OnePerCore>("Time: %f Value: %f", timeTemp, valueTemp);
          timeValuePairs[timeTemp] = valueTemp;
        }

        datafile.close();
        // the default iterator for maps traverses in key order, so no sort is needed.

        std::vector<PhysicalTime> times(0);
        std::vector<PhysicalSpeed> values(0);

        // Must convert into vectors since LinearInterpolate works on a pair of vectors
        PhysicalTime t_max = endTS * timeStepLength;
        for (auto& [time, speed]: timeValuePairs)
        {
            /* If the time value in the input file stretches BEYOND the end of the simulation, then insert an interpolated end value and exit the loop. */
            if (time > t_max) {
                PhysicalTime time_diff = t_max - times.back();

                PhysicalTime time_diff_ratio = time_diff / (time - times.back());
                PhysicalSpeed vel_diff = speed - values.back();

                PhysicalSpeed final_velocity = values.back() + time_diff_ratio * vel_diff;

                times.push_back(t_max);
                values.push_back(final_velocity);
                break;
            }

            times.push_back(time);
            values.push_back(speed);
        }

        /* If the time values in the input file end BEFORE the planned end of the simulation, then loop the profile afterwards (using %TimeStepsInInletVelocityProfile). */
        PhysicalTime duration = (times.back() - times.front());
        int TimeStepsInInletVelocityProfile = duration / timeStepLength;

        // Check if last point's value matches the first
        if (values.back() != values.front())
          throw Exception() << "Last point's value does not match the first point's value in "
              << velocityFilePath;

        // extend the table to one past the total time steps, so that the table is valid in the end-state, where the zero indexed time step is equal to the limit.
        //auto totalTimeSteps = endTS - startTS;
        velocityTable.resize(endTS + 1);
        // Now convert these vectors into arrays using linear interpolation
        for (unsigned int time_step = 0; time_step <= endTS; ++time_step)
        {
            // The "% TimeStepsInInletVelocityProfile" here is to prevent profile stretching (it will loop instead).
            double point = times.front() +
                    (double(time_step % TimeStepsInInletVelocityProfile) / double(endTS)) * duration;

            PhysicalSpeed vel = LinearInterpolate(times, values, point);
            velocityTable[time_step] = units->ConvertVelocityToLatticeUnits(vel);
        }

      }

    LatticeVelocity InOutLetFileVelocity::GetVelocity(const LatticePosition& x,
                                                      const LatticeTimeStep t) const
    {
        if (!useWeightsFromFile)
        {
          // v(r) = vMax (1 - r**2 / a**2)
          // where r is the distance from the centreline
          LatticePosition displ = x - position;
          LatticeDistance z = Dot(displ, normal);
          Dimensionless rSqOverASq = (displ.GetMagnitudeSquared() - z * z) / (radius * radius);
          HASSERT(rSqOverASq <= 1.0);

          // Get the max velocity
          LatticeSpeed max = velocityTable[t];

          // Brackets to ensure that the scalar multiplies are done before vector * scalar.
          return normal * (max * (1. - rSqOverASq));
        }
        else
        {
          /* These absolute normal values can still be negative here,
           * but are corrected below to become positive.
           * */
          auto abs_normal = normal;

          /* Prevent division by 0 errors if the normals are 0.0. */
          for (auto& comp: abs_normal) {
              comp = std::max(comp, 0.0000001);
          }

          /*bool logging = false;
          if (402.9 < x.x && x.x < 403.1 && 312.9 < x.y && x.y < 313.1 && 160.4 < x.z && x.z < 160.6)
           {
           logging = true;
           }

          if (logging)
          {
            log::Logger::Log<log::Warning, log::OnePerCore>("%f %f %f", x.x, x.y, x.z);
          }*/

          int xyz_directions[3] = { 1, 1, 1 };

          std::vector<int> xyz;
          xyz.push_back(0);
          xyz.push_back(0);
          xyz.push_back(0);

          double xyz_residual[3] = {0.0, 0.0, 0.0};
          /* The residual values increase by the normal values at every time step. When they hit >1.0, then
           * xyz is incremented and a new grid point is attempted.
           * In addition, the specific residual value is decreased by 1.0. */

          for (int i = 0; i < 3; ++i)  {
              if (normal[i] < 0.0) {
                  xyz_directions[i] = -1;
                  xyz[i] = floor(x[i]);
                  abs_normal[i] = -abs_normal[i];
                  /* Start with a negative residual because we already moved partially in this direction. */
                  xyz_residual[i] = -(x[i] - floor(x[i]));
              } else {
                  xyz[i] = std::ceil(x[i]);
                  xyz_residual[i] = -(std::ceil(x[i]) - x[i]);
              }
          }

	  auto v_tot = LatticeVelocity::Zero();
          int iterations = 0;

          while (iterations < 3)
          {
            if (weights_table.count(xyz) > 0)
            {
              v_tot = normal * weights_table.at(xyz) * velocityTable[t];
              //log::Logger::Log<log::Warning, log::OnePerCore>("%f %f %f %f",
              //                                                              x.x,
              //                                                              x.y,
              //                                                              x.z, v_tot);
              return v_tot;
            }

            /*if (logging)
            {
              log::Logger::Log<log::Warning, log::OnePerCore>("%f %f %f %d %d %d",
                                                              x.x,
                                                              x.y,
                                                              x.z,
                                                              xyz[0],
                                                              xyz[1],
                                                              xyz[2]);
            }*/

            /* Propagate residuals to the move to the next grid point. */
            double xstep = (1.0 - xyz_residual[0]) / abs_normal[0];
            double ystep = (1.0 - xyz_residual[1]) / abs_normal[1];
            double zstep = (1.0 - xyz_residual[2]) / abs_normal[2];

            //log::Logger::Log<log::Warning, log::OnePerCore>("%f %f %f", xstep, ystep, zstep);

            double all_step = 0.0;
            int xyz_change = 0;

            if(xstep < ystep) {
              if (xstep < zstep) {
                all_step = xstep;
                xyz_change = 0;
              } else {
                if (ystep < zstep) {
                  all_step = ystep;
                  xyz_change = 1;
                } else {
                  all_step = zstep;
                  xyz_change = 2;
                }
              }

            } else {
              if (ystep < zstep) {
                all_step = ystep;
                xyz_change = 1;
              } else {
                all_step = zstep;
                xyz_change = 2;
              }
            }

            xyz_residual[0] += abs_normal[0] * all_step;
            xyz_residual[1] += abs_normal[1] * all_step;
            xyz_residual[2] += abs_normal[2] * all_step;

            xyz[xyz_change] += xyz_directions[xyz_change];

            //if(xyz_residual[xyz_change] < 1.0) {
            //  log::Logger::Log<log::Error, log::Singleton>("ERROR: Residual bug in vInlet: %f %f %f %f", x.x, x.y, x.z, xyz_residual[xyz_change]);
            //}

            xyz_residual[xyz_change] -= 1.0;

            iterations++;
          }

          /* Lists the sites which should be in the wall, outside of the main inlet.
           * If you are unsure, you can increase the log level of this, run HemeLb
           * for 1 time step, and plot these points out. */
          log::Logger::Log<log::Trace, log::OnePerCore>("%f %f %f", x.x(), x.y(), x.z());
          return normal * 0.0;
        }

      }

      void InOutLetFileVelocity::Initialise(const util::UnitConverter* unitConverter)
      {
        log::Logger::Log<log::Warning, log::OnePerCore>("Initializing vInlet.");
        units = unitConverter;

        if constexpr (build_info::USE_VELOCITY_WEIGHTS_FILE) {
            useWeightsFromFile = true;
        } else {
            useWeightsFromFile = false;
        }

        if(useWeightsFromFile) {
          //if the new velocity approximation is enabled, then we want to create a lookup table here.
          const std::string in_name = velocityFilePath + ".weights.txt";
          if (!std::filesystem::exists(in_name))
              throw Exception() << "File does not exist: " << in_name;

          /* Load and read file. */
          std::fstream myfile;
          myfile.open(in_name, std::ios_base::in);
          log::Logger::Log<log::Warning, log::OnePerCore>("Loading weights file: %s",
                                                        in_name.c_str());

          std::string input_line;
          /* input files are in ASCII, in format:
           *
           * coord_x coord_y coord_z weights_value
           *
           * */
          while (myfile.good()) //(std::getline(myfile, input_line))
          {
            int x, y, z;
            double v;
            myfile >> x >> y >> z >> v;

            std::vector<int> xyz;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
            weights_table[xyz] = v;

            log::Logger::Log<log::Trace, log::OnePerCore>("%lld %lld %lld %f",
            x,
            y,
            z,
            weights_table[xyz]);
          }
          myfile.close();
        }
      }

}
