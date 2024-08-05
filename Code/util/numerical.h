// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_NUMERICAL_H
#define HEMELB_UTIL_NUMERICAL_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

namespace hemelb::util
{

    /**
     * Helper for performing x^n when n is an integer. (And ought to be more efficient than
     * using pow in cmath).
     *
     * Note that this function doesn't check for bad maths like 0^0.
     *
     * @param x
     * @param n
     * @return
     */
    template<typename T>
    constexpr T IntegerPower(T x, long n)
    {
        if (n == 0) {
            return T(1);
        } else if (n > 0) {
            return x * IntegerPower(x, n - 1);
        } else {
            return IntegerPower(x, n + 1) / x;
        }
    }



    template<typename T>
    constexpr bool IsInRange(T x, T min, T max)
    {
        return ( (x >= min) && (x <= max));
    }

    /**
     * The Newton_Raphson method takes in a functor with the operator () overloaded, which
     * should return void and take in three doubles: x, f, df. x is the variable to
     * be solved for, f and df are the function and derivative values at x respectively.
     * The function calculates f and df at the given x and stores in the given f and df.
     * The other two arguments for NewtonRaphson are the initial guess and desired accuracy.
     */
    template<class F>
    double NewtonRaphson(F&& func, double x0, double alphaAcc)
    {
        double x = x0, dx;
        double f, df;

        for (int i = 0; i < 20; i++) {
            func(x, f, df);

            dx = f / df;
            x -= dx;

            if (std::fabs(dx) < alphaAcc)
            {
                return x;
            }
        }

        /*
         * TODO: Implement some sensible way of dealing with too many iterations
         */

        return x;
    }

    /*
     * Finds root using Brent's method. Needs to be given a bracket enclosing the root.
     * The caller must check if a root is enclosed so that he can specify the result in that case
     * Since it must check for this it will have the values of the function at those points
     * so they need to be passed on as well.
     */
    template<class F>
    double Brent(F&& func, double xLowerIn, double yLowerIn, double xHigherIn,
                 double yHigherIn, double xAccuracy, double yAccuracy)
    {
        double xLower = xLowerIn, yLower = yLowerIn;
        double xHigher = xHigherIn, yHigher = yHigherIn;
        double xBoundNew, yBoundNew;
        double xBoundOld; // First set after first iteration hence mflag
        double xSolution = xHigher, ySolution = yHigher;

        if (std::fabs(yLower) < std::fabs(yHigher))
        {
            double temp = yLower;
            yLower = yHigher;
            yHigher = temp;
            temp = xLower;
            xLower = xHigher;
            xHigher = temp;
        }

        xBoundNew = xLower;
        yBoundNew = yLower;

        bool mflag = true;

        while (std::fabs(xHigher - xLower) > xAccuracy && std::fabs(yHigher) > yAccuracy
               && std::fabs(ySolution) > yAccuracy)
        {
            if (yLower != yBoundNew && yHigher != yBoundNew)
            {
                xSolution = (xLower * yHigher * yBoundNew)
                            / ( (yLower - yHigher) * (yLower - yBoundNew))
                            + (xHigher * yLower * yBoundNew) / ( (yHigher - yLower) * (yHigher - yBoundNew))
                            + (xBoundNew * yLower * yHigher)
                              / ( (yBoundNew - yLower) * (yBoundNew - yHigher));
            }
            else
            {
                xSolution = xHigher - yHigher * (xHigher - xLower) / (yHigher - yLower);
            }

            // s is not between (3a + b)/4 and b
            bool condition1 = (xLower < xHigher ?
                               (xSolution < (3 * xLower + xHigher) / 4.0 || xSolution > xHigher) :
                               (xSolution > (3 * xLower + xHigher) / 4.0 || xSolution < xHigher)); // mflag is set and |s−b| ≥ |b−c| / 2)
            bool condition2 = mflag
                              && std::fabs(xSolution - xHigher) >= std::fabs(xHigher - xBoundNew) / 2.0;
            // mflag is cleared and |s−b| ≥ |c−d| / 2
            bool condition3 = !mflag
                              && std::fabs(xSolution - xHigher) >= std::fabs(xBoundNew - xBoundOld) / 2.0;
            // mflag is set and |b−c| < |δ|
            bool condition4 = mflag && std::fabs(xHigher - xBoundNew) < xAccuracy;
            // mflag is cleared and |c−d| < |δ|
            bool condition5 = !mflag && std::fabs(xBoundNew - xBoundOld) < xAccuracy;

            if (condition1 || condition2 || condition3 || condition4 || condition5)
            {
                xSolution = (xLower + xHigher) / 2.0;
                mflag = true;
            }
            else
            {
                mflag = false;
            }

            func(xSolution, ySolution);
            xBoundOld = xBoundNew;
            xBoundNew = xHigher;
            yBoundNew = yHigher;

            if (yLower * ySolution < 0)
            {
                xHigher = xSolution;
                yHigher = ySolution;
            }
            else
            {
                xLower = xSolution;
                yLower = ySolution;
            }

            if (std::fabs(yLower) < std::fabs(yHigher))
            {
                double temp = yLower;
                yLower = yHigher;
                yHigher = temp;
                temp = xLower;
                xLower = xHigher;
                xHigher = temp;
            }
        }

        if (std::fabs(yHigher) < std::fabs(ySolution))
            return xHigher;
        else
            return xSolution;
    }

}

#endif
