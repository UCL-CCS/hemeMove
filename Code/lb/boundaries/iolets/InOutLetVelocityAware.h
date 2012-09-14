/*
 * InOutLetVelocityAware.h
 *
 *  Created on: 26 Jun 2012
 *      Author: derek
 */

#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETVELOCITYAWARE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETVELOCITYAWARE_H

#include "lb/boundaries/iolets/InOutLetMultiscale.h"
#include "geometry/LatticeData.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "log/Logger.h"
//#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    class MacroscopicPropertyCache;

    namespace boundaries
    {
      namespace iolets
      {

        /***
         * An inlet/outlet whose density is obtained from the multiscale intercommunicator
         * We envisage communication of velocity information outwards to other processes
         * the velocity SharedValue is a place-holder for this.
         * We set 0.1 as an arbitrary value for the output velocity for now.
         * We do not yet have an understanding of the necessary physics, or the necessary computational infrastructure,
         * for building boundary conditions which have access to their nearby velocity field.
         * The min and max pressure SharedValues are placeholders for information needed by the steering and visualisation code.
         * The steering and visualisation code requires minimum and maximum pressure values.
         */
        class InOutLetVelocityAware : public InOutLetMultiscale
        {
          public:
            InOutLetVelocityAware() :
                InOutLetMultiscale(), NDM(NULL), propertyCache(NULL), sitesWhichNeighbourThisBoundary()

            {
            }
            /***ss
             * The shared values are registered through the initialiser-list syntactic sugar.
             */
            InOutLetVelocityAware(const InOutLetVelocityAware &other) :
                InOutLetMultiscale(other), NDM(other.NDM), propertyCache(other.propertyCache), sitesWhichNeighbourThisBoundary()
            {
              /* Add a velocity aware exchange here if needed?*/
            }

            void InitialiseNeighbouringSites(geometry::neighbouring::NeighbouringDataManager *manager,
                                             geometry::LatticeData * latDat,
                                             hemelb::lb::MacroscopicPropertyCache* globalPropertyCache,
                                             std::vector<site_t> invBList);

            PhysicalVelocity GetVelocity() const;
            virtual ~InOutLetVelocityAware()
            {
            }

            //virtual void DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig);

            virtual InOutLet* Clone() const
            {
              InOutLetVelocityAware* copy = new InOutLetVelocityAware(*this);
              return copy;
            }
          protected:
            geometry::neighbouring::NeighbouringDataManager *NDM;
            std::vector<site_t> sitesWhichNeighbourThisBoundary;
            geometry::LatticeData * latticeData;
            hemelb::lb::MacroscopicPropertyCache* propertyCache;
        }
        ;
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETVELOCITYAWARE_H */
