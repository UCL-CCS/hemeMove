#ifndef HEMELB_UNITTESTS_EXTRACTION_GEOMETRYSELECTORTESTS_H
#define HEMELB_UNITTESTS_EXTRACTION_GEOMETRYSELECTORTESTS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <cppunit/TestFixture.h>
#include "extraction/LbDataSourceIterator.h"
#include "extraction/StraightLineGeometrySelector.h"
#include "extraction/PlaneGeometrySelector.h"
#include "extraction/WholeGeometrySelector.h"
#include "unittests/FourCubeLatticeData.h"

namespace hemelb
{
  namespace unittests
  {
    namespace extraction
    {
      class GeometrySelectorTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(GeometrySelectorTests);
          CPPUNIT_TEST(TestStraightLineGeometrySelector);
          CPPUNIT_TEST(TestPlaneGeometrySelector);
          CPPUNIT_TEST(TestWholeGeometrySelector);CPPUNIT_TEST_SUITE_END();

        public:
          GeometrySelectorTests() :
              VoxelSize(0.01), CubeSize(10), CentreCoordinate( ((distribn_t) CubeSize - 1.0) / 2.0), planeNormal(1.0),

              planePosition(CentreCoordinate * VoxelSize), planeRadius(distribn_t(CubeSize)
                  * VoxelSize / 3.0), lineEndPoint1(CentreCoordinate * VoxelSize), lineEndPoint2( (CubeSize
                  + 1) * VoxelSize)
          {

          }

          void setUp()
          {
            latticeData = FourCubeLatticeData::Create(CubeSize, 1);
            simState = new lb::SimulationState(1000, 5);
            propertyCache = new hemelb::lb::MacroscopicPropertyCache(*simState, *latticeData);
            lbmParams = new lb::LbmParameters(0.1, VoxelSize);
            unitConverter = new hemelb::util::UnitConverter(lbmParams, simState, VoxelSize);
            dataSourceIterator = new hemelb::extraction::LbDataSourceIterator(*propertyCache,
                                                                              *latticeData,
                                                                              *unitConverter);

            planeGeometrySelector = new hemelb::extraction::PlaneGeometrySelector(planePosition,
                                                                                  planeNormal);
            planeGeometrySelectorWithRadius =
                new hemelb::extraction::PlaneGeometrySelector(planePosition,
                                                              planeNormal,
                                                              planeRadius);
            straightLineGeometrySelector = new hemelb::extraction::StraightLineGeometrySelector(lineEndPoint1,
                                                                                lineEndPoint2);
            wholeGeometrySelector = new hemelb::extraction::WholeGeometrySelector();
          }

          void tearDown()
          {
            delete planeGeometrySelector;
            delete planeGeometrySelectorWithRadius;
            delete straightLineGeometrySelector;
            delete wholeGeometrySelector;

            delete dataSourceIterator;
            delete unitConverter;
            delete lbmParams;
            delete propertyCache;
            delete simState;
            delete latticeData;
          }

          void TestStraightLineGeometrySelector()
          {
            TestOutOfGeometrySites(straightLineGeometrySelector);

            // The line runs from the centre to the (size,size,size) point.
            // This includes anything with all three coordinates the same, above
            // a minimum.
            const site_t coordMin = (CubeSize) / 2;

            // Gather the list of coordinates we expect to be on the line.
            std::vector<util::Vector3D<site_t> > includedCoords;

            for (site_t xCoord = coordMin; xCoord < CubeSize; ++xCoord)
            {
              includedCoords.push_back(util::Vector3D<site_t>(xCoord));
            }

            TestExpectedIncludedSites(straightLineGeometrySelector, includedCoords);
          }

          void TestPlaneGeometrySelector()
          {
            TestOutOfGeometrySites(planeGeometrySelector);
            TestOutOfGeometrySites(planeGeometrySelectorWithRadius);

            // Gather the list of coordinates we expect to be on the plane with and
            // without a radius used.
            std::vector<util::Vector3D<site_t> > includedCoordsWithRadius;
            std::vector<util::Vector3D<site_t> > includedCoordsWithoutRadius;

            // The plane has normal (1, 1, 1), is centred in the cube and has a radius of
            // CubeSize / 3 lattice units (when the radius is in use).
            for (site_t xCoord = 0; xCoord < CubeSize; ++xCoord)
            {
              for (site_t yCoord = 0; yCoord < CubeSize; ++yCoord)
              {
                for (site_t zCoord = 0; zCoord < CubeSize; ++zCoord)
                {
                  // Use that p.n = x.n for x on the same plane.
                  // I.e. the current point's coordinate dotted with the normal must be roughly equal
                  // to the centre point's coordinate dotted with the normal for the current point
                  // to be on the plane.
                  if (std::abs(distribn_t(xCoord + yCoord + zCoord) - 3 * CentreCoordinate) > 0.5)
                  {
                    continue;
                  }

                  includedCoordsWithoutRadius.push_back(util::Vector3D<site_t>(xCoord,
                                                                               yCoord,
                                                                               zCoord));

                  // Compute the distance from the centre point, include the site if it is within the radius.
                  if ( (util::Vector3D<distribn_t>(xCoord, yCoord, zCoord)
                      - util::Vector3D<distribn_t>(CentreCoordinate)).GetMagnitude()
                      < distribn_t(CubeSize) / 3.0)
                  {
                    includedCoordsWithRadius.push_back(util::Vector3D<site_t>(xCoord,
                                                                              yCoord,
                                                                              zCoord));
                  }
                }
              }
            }

            TestExpectedIncludedSites(planeGeometrySelector, includedCoordsWithoutRadius);
            TestExpectedIncludedSites(planeGeometrySelectorWithRadius, includedCoordsWithRadius);
          }

          void TestWholeGeometrySelector()
          {
            TestOutOfGeometrySites(wholeGeometrySelector);

            // Variables to count the number of sites encountered, and store the iterated-over
            // positions, and properties.
            int count = 0;
            util::Vector3D<site_t> position;
            util::Vector3D<float> velocity;
            float pressure, stress;

            dataSourceIterator->Reset();
            while (dataSourceIterator->ReadNext(position, pressure, velocity, stress))
            {
              ++count;
              // Every site returned should be included.
              CPPUNIT_ASSERT(wholeGeometrySelector->Include(*dataSourceIterator, position));
            }

            // The number of sites passed by the iterator should be the whole cube.
            CPPUNIT_ASSERT_EQUAL(CubeSize * CubeSize * CubeSize, count);
          }

        private:
          void TestOutOfGeometrySites(hemelb::extraction::GeometrySelector* geometrySelector)
          {
            util::Vector3D<site_t> invalidLocation(-1);
            CPPUNIT_ASSERT(!geometrySelector->Include(*dataSourceIterator, invalidLocation));

            invalidLocation = util::Vector3D<site_t>(CubeSize);
            CPPUNIT_ASSERT(!geometrySelector->Include(*dataSourceIterator, invalidLocation));
          }

          void TestExpectedIncludedSites(hemelb::extraction::GeometrySelector* geometrySelector,
                                         std::vector<util::Vector3D<site_t> >& includedSites)
          {
            // Variables to store the iterated-over positions, and properties.
            util::Vector3D<site_t> position;
            util::Vector3D<float> velocity;
            float pressure, stress;

            dataSourceIterator->Reset();
            while (dataSourceIterator->ReadNext(position, pressure, velocity, stress))
            {
              bool expectedIncluded = std::count(includedSites.begin(),
                                                 includedSites.end(),
                                                 position) > 0;

              std::stringstream msg;

              msg << "Site at " << position.x << "," << position.y << "," << position.z << " was ";

              if (!expectedIncluded)
              {
                msg << "not ";
              }

              msg << "expected to be included but actually was";

              if (!geometrySelector->Include(*dataSourceIterator, position))
              {
                msg << " not.";
              }
              else
              {
                msg << ".";
              }

              CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(),
                                           expectedIncluded,
                                           geometrySelector->Include(*dataSourceIterator,
                                                                     position));
            }
          }

          const double VoxelSize;
          const int CubeSize;
          const distribn_t CentreCoordinate;

          const util::Vector3D<distribn_t> planeNormal;
          const util::Vector3D<distribn_t> planePosition;
          const distribn_t planeRadius;
          const util::Vector3D<distribn_t> lineEndPoint1;
          const util::Vector3D<distribn_t> lineEndPoint2;

          unittests::FourCubeLatticeData* latticeData;
          lb::SimulationState* simState;
          hemelb::lb::MacroscopicPropertyCache* propertyCache;
          hemelb::lb::LbmParameters* lbmParams;
          hemelb::util::UnitConverter* unitConverter;
          hemelb::extraction::LbDataSourceIterator* dataSourceIterator;

          hemelb::extraction::PlaneGeometrySelector* planeGeometrySelector;
          hemelb::extraction::PlaneGeometrySelector* planeGeometrySelectorWithRadius;
          hemelb::extraction::StraightLineGeometrySelector* straightLineGeometrySelector;
          hemelb::extraction::WholeGeometrySelector* wholeGeometrySelector;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(GeometrySelectorTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_EXTRACTION_GEOMETRYSELECTORTESTS_H */
