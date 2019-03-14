
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMCONFIG_H
#define HEMELB_CONFIGURATION_SIMCONFIG_H

#include <vector>
#include "util/Vector3D.h"
#include "lb/LbmParameters.h"
#include "lb/iolets/InOutLets.h"
#include "lb/stents/Stents.h"
#include "extraction/PropertyOutputFile.h"
#include "extraction/GeometrySelectors.h"
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
  namespace configuration
  {
    template<typename T>
    void GetDimensionalValue(const io::xml::Element& elem, const std::string& units, T& value)
    {
      const std::string& got = elem.GetAttributeOrThrow("units");
      if (got != units)
      {
        throw Exception() << "Invalid units for element " << elem.GetPath() << ". Expected '"
            << units << "', got '" << got << "'";
      }

      elem.GetAttributeOrThrow("value", value);
    }

    class SimConfig
    {
      public:
        /**
         * Bundles together various configuration parameters concerning simulation monitoring
         */
        struct MonitoringConfig
        {
            MonitoringConfig() :
                doConvergenceCheck(false), convergenceRelativeTolerance(0), convergenceTerminate(false),
                    doIncompressibilityCheck(false)
            {
            }
            bool doConvergenceCheck; ///< Whether to turn on the convergence check or not
            extraction::OutputField::FieldType convergenceVariable; ///< Macroscopic variable used to check for convergence
            double convergenceReferenceValue; ///< Reference value used to normalise an absolute error (making it relative)
            double convergenceRelativeTolerance; ///< Convergence check relative tolerance
            bool convergenceTerminate; ///< Whether to terminate a converged run or not
            bool doIncompressibilityCheck; ///< Whether to turn on the IncompressibilityChecker or not
        };

        static SimConfig* New(const std::string& path);

      protected:
        SimConfig(const std::string& path);
        void Init();

      public:
        virtual ~SimConfig();

        void Save(std::string path); // TODO this method should be able to be CONST
        // but because it uses DoIo, which uses one function signature for both reading and writing, it cannot be.

        const util::Vector3D<float> & GetVisualisationCentre() const
        {
          return visualisationCentre;
        }
        const std::vector<lb::iolets::InOutLet*> & GetInlets() const
        {
          return inlets;
        }
        const std::vector<lb::iolets::InOutLet*> & GetOutlets() const
        {
          return outlets;
        }
        const std::vector<lb::stents::Stent*> & GetStents() const
        {
          return stents;
        }
        lb::StressTypes GetStressType() const
        {
          return stressType;
        }
        float GetVisualisationLongitude() const
        {
          return visualisationLongitude;
        }
        float GetVisualisationLatitude() const
        {
          return visualisationLatitude;
        }
        float GetVisualisationZoom() const
        {
          return visualisationZoom;
        }
        float GetVisualisationBrightness() const
        {
          return visualisationBrightness;
        }
        float GetMaximumVelocity() const
        {
          return maxVelocity;
        }
        float GetMaximumStress() const
        {
          return maxStress;
        }
        const std::string & GetDataFilePath() const
        {
          return dataFilePath;
        }
        LatticeTimeStep GetTotalTimeSteps() const
        {
          return totalTimeSteps;
        }
        LatticeTimeStep GetWarmUpSteps() const
        {
          return warmUpSteps;
        }
        PhysicalTime GetTimeStepLength() const
        {
          return timeStepSeconds;
        }
        PhysicalDistance GetVoxelSize() const
        {
          return voxelSizeMetres;
        }
        PhysicalPosition GetGeometryOrigin() const
        {
          return geometryOriginMetres;
        }
        unsigned int PropertyOutputCount() const
        {
          return propertyOutputs.size();
        }
        extraction::PropertyOutputFile * GetPropertyOutput(unsigned int index) const
        {
          return propertyOutputs[index];
        }
        std::vector<extraction::PropertyOutputFile*> const GetPropertyOutputs() const
        {
          return propertyOutputs;
        }
        const std::string GetColloidConfigPath() const
        {
          return colloidConfigPath;
        }
        /**
         * True if the XML file has a section specifying colloids.
         * @return
         */
        bool HasColloidSection() const;

        /**
         * Returns the pressure to be used to initialise all the fluid sites in the domain
         * @return initial pressure
         */
        LatticeDensity GetInitialPressure() const;

        LatticeDensity GetInitialDensity() const;

        const util::UnitConverter& GetUnitConverter() const;

        /**
         * Return the configuration of various checks/test
         * @return monitoring configuration
         */
        const MonitoringConfig* GetMonitoringConfiguration() const;

      protected:
        /**
         * Create the unit converter - virtual so that mocks can override it.
         */
        virtual void CreateUnitConverter();

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                            const std::string& requiredBC);

        template<typename T>
        void GetDimensionalValueInLatticeUnits(const io::xml::Element& elem,
                                               const std::string& units, T& value)
        {
          GetDimensionalValue(elem, units, value);
          value = unitConverter->ConvertToLatticeUnits(units, value);
        }

        template<typename T>
        T GetDimensionalValueInLatticeUnits(const io::xml::Element& elem, const std::string& units)
        {
          T ans;
          GetDimensionalValueInLatticeUnits(elem, units, ans);
          return ans;
        }

      private:
        void DoIO(const io::xml::Element xmlNode);
        void DoIOForSimulation(const io::xml::Element simEl);
        void DoIOForGeometry(const io::xml::Element geometryEl);

        std::vector<lb::iolets::InOutLet*> DoIOForInOutlets(const io::xml::Element xmlNode);

        std::vector<hemelb::lb::stents::Stent*> DoIOForStents(const io::xml::Element xmlNode);

        void DoIOForBaseInOutlet(const io::xml::Element& ioletEl, lb::iolets::InOutLet* value);

        lb::iolets::InOutLet* DoIOForPressureInOutlet(const io::xml::Element& ioletEl);
        lb::stents::Stent* DoIOForDensityStent(const io::xml::Element& stentEl);
        lb::iolets::InOutLetCosine* DoIOForCosinePressureInOutlet(const io::xml::Element& ioletEl);
        lb::stents::StentConstant* DoIOForConstantDensityStent(const io::xml::Element& stentEl);
        lb::iolets::InOutLetFile* DoIOForFilePressureInOutlet(const io::xml::Element& ioletEl);
        lb::stents::StentFile* DoIOForFileDensityStent(const io::xml::Element& stentEl);
        lb::iolets::InOutLetMultiscale* DoIOForMultiscalePressureInOutlet(
            const io::xml::Element& ioletEl);

        lb::iolets::InOutLet* DoIOForVelocityInOutlet(const io::xml::Element& ioletEl);
        lb::stents::Stent* DoIOForFluxStent(const io::xml::Element& stentEl);
        lb::iolets::InOutLetParabolicVelocity* DoIOForParabolicVelocityInOutlet(
            const io::xml::Element& ioletEl);
        lb::iolets::InOutLetHollowParabolicVelocity* DoIOForHollowParabolicVelocityInOutlet(
            const io::xml::Element& ioletEl);
        lb::stents::StentConstantFlux* DoIOForConstantFluxStent(
            const io::xml::Element& stentEl);
        /**
         * Reads a Womersley velocity iolet definition from the XML config file and returns
         * an InOutLetWomersleyVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetWomersleyVelocity object
         */
        lb::iolets::InOutLetWomersleyVelocity* DoIOForWomersleyVelocityInOutlet(
            const io::xml::Element& ioletEl);

        /**
         * Reads a file velocity iolet definition from the XML config file and returns
         * an InOutLetFileVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetFileVelocity object
         */
        lb::iolets::InOutLetFileVelocity* DoIOForFileVelocityInOutlet(
            const io::xml::Element& ioletEl);
        lb::iolets::InOutLetHollowFileVelocity* DoIOForHollowFileVelocityInOutlet(
            const io::xml::Element& ioletEl);
        lb::stents::StentFileFlux* DoIOForFileFluxStent(
            const io::xml::Element& stentEl);

        lb::stents::StentCoating* DoIOForCoatingStent(const io::xml::Element& stentEl);

        void DoIOForProperties(const io::xml::Element& xmlNode);
        void DoIOForProperty(io::xml::Element xmlNode, bool isLoading);
        extraction::OutputField DoIOForPropertyField(const io::xml::Element& xmlNode);
        extraction::PropertyOutputFile* DoIOForPropertyOutputFile(
            const io::xml::Element& propertyoutputEl);
        extraction::StraightLineGeometrySelector* DoIOForLineGeometry(
            const io::xml::Element& xmlNode);
        extraction::PlaneGeometrySelector* DoIOForPlaneGeometry(const io::xml::Element&);
        extraction::SurfacePointSelector* DoIOForSurfacePoint(const io::xml::Element&);

        void DoIOForInitialConditions(io::xml::Element parent);
        void DoIOForVisualisation(const io::xml::Element& visEl);

        /**
         * Reads monitoring configuration from XML file
         *
         * @param monEl in memory representation of <monitoring> xml element
         */
        void DoIOForMonitoring(const io::xml::Element& monEl);

        /**
         * Reads configuration of steady state flow convergence check from XML file
         *
         * @param convEl in memory representation of the <steady_flow_convergence> XML element
         */
        void DoIOForSteadyFlowConvergence(const io::xml::Element& convEl);

        /**
         * Reads the configuration of one of the possible several converge criteria provided
         *
         * @param criterionEl in memory representation of the <criterion> XML element
         */
        void DoIOForConvergenceCriterion(const io::xml::Element& criterionEl);

        const std::string& xmlFilePath;
        io::xml::Document* rawXmlDoc;
        std::string dataFilePath;

        util::Vector3D<float> visualisationCentre;
        float visualisationLongitude;
        float visualisationLatitude;
        float visualisationZoom;
        float visualisationBrightness;
        float maxVelocity;
        float maxStress;
        lb::StressTypes stressType;
        std::vector<extraction::PropertyOutputFile*> propertyOutputs;
        std::string colloidConfigPath;
        /**
         * True if the file has a colloids section.
         */
        bool hasColloidSection;
        PhysicalPressure initialPressure_mmHg; ///< Pressure used to initialise the domain
        PhysicalDensity initialDensity_Kg_per_m3;
        MonitoringConfig monitoringConfig; ///< Configuration of various checks/tests

      protected:
        // These have to contain pointers because there are multiple derived types that might be
        // instantiated.
        std::vector<lb::iolets::InOutLet*> inlets;
        std::vector<lb::iolets::InOutLet*> outlets;
        std::vector<lb::stents::Stent*> stents;
        PhysicalTime timeStepSeconds;
        unsigned long totalTimeSteps;
        unsigned long warmUpSteps;
        PhysicalDistance voxelSizeMetres;
        PhysicalPosition geometryOriginMetres;
        util::UnitConverter* unitConverter;
    };
  }
}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
