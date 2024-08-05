// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/iolets/InOutLets.h"
#include "configuration/SimConfigReader.h"
#include "configuration/SimBuilder.h"
#include "resources/Resource.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb::tests
{
    //using namespace lb::iolets;
    using resources::Resource;

    class UncheckedSimConfigReader : public configuration::SimConfigReader {
    public:
        using configuration::SimConfigReader::SimConfigReader;
    protected:
        void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                            std::string_view requiredBC) const override
        {

        }
    };

    class ConcreteIolet : public lb::InOutLet {
        [[nodiscard]] InOutLet* clone() const override
        {
            return new ConcreteIolet(*this);
        }
        LatticeDensity GetDensityMin() const override
        {
            return 1.0;
        }
        LatticeDensity GetDensityMax() const override
        {
            return 1.0;
        }
        LatticeDensity GetDensity(hemelb::LatticeTimeStep) const override
        {
            return 1.0;
        }
        void Reset(hemelb::lb::SimulationState&) override
        {
        }
    };

    TEST_CASE_METHOD(helpers::FolderTestFixture, "InOutLetTests") {

        SECTION("TestCosineConstruct") {
            // Bootstrap ourselves an iolet, by loading config.xml.
            UncheckedSimConfigReader reader(Resource("config.xml").Path());
            auto config = reader.Read();
            configuration::SimBuilder builder(config);
            auto& converter = *builder.GetUnitConverter();

            auto cosineConfig = config.GetInlets()[0];
            auto cosine = util::clone_dynamic_cast<lb::InOutLetCosine>(builder.BuildIolet(cosineConfig));

            lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                            config.GetTotalTimeSteps());
            // at this stage, Initialise() has not been called, so the
            // unit converter will be invalid, so we will not be able to
            // convert to physical units.
            cosine->Initialise(&converter);
            cosine->Reset(state);

            // Check the cosine IOLET contains the values expected given the file.
            /*
             * <inlet>
             <pressure amplitude="0.0" mean="80.1" phase="0.0" period="0.6"/>
             <normal x="0.0" y="0.0" z="1.0" />
             <position x="-1.66017717834e-05" y="-4.58437586355e-05" z="-0.05" />
             </inlet>
            */
            REQUIRE(Approx(80.1) == converter.ConvertPressureToPhysicalUnits(cosine->GetPressureMean()));
            REQUIRE(0.0 == cosine->GetPressureAmp());
            REQUIRE(0.0 == cosine->GetPhase());
            REQUIRE(Approx(0.6) == converter.ConvertTimeToPhysicalUnits(cosine->GetPeriod()));
            auto expected = ApproxVector<PhysicalPosition>{-1.66017717834e-05, -4.58437586355e-05, -0.05};
            PhysicalPosition actual =
                    converter.ConvertPositionToPhysicalUnits(cosine->GetPosition());
            REQUIRE(expected == actual);

            REQUIRE(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0) == cosine->GetNormal());

            // Set an appropriate target value for the density, the maximum.
            LatticeDensity targetMeanDensity = converter.ConvertPressureToLatticeUnits(80.1) / Cs2;
            // Check that the cosine formula correctly produces mean value
            REQUIRE(Approx(targetMeanDensity) == cosine->GetDensityMean());
            REQUIRE(Approx(targetMeanDensity) == cosine->GetDensity(0));
        }

        SECTION("TestFileConstruct") {
            // Bootstrap ourselves a file inlet, by loading an appropriate
            // config file.  We have to move to a tempdir, as the path
            // from the inlet to the iolet.txt file is a relative path
            CopyResourceToTempdir("iolet.txt");
            MoveToTempdir();

            UncheckedSimConfigReader reader(Resource("config_file_inlet.xml").Path());
            auto config = reader.Read();
            configuration::SimBuilder builder(config);
            auto& converter = *builder.GetUnitConverter();

            lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                            config.GetTotalTimeSteps());
            auto fileConfig = config.GetInlets()[0];
            auto file = util::clone_dynamic_cast<lb::InOutLetFile>(builder.BuildIolet(fileConfig));
            // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
            file->Initialise(&converter);
            file->Reset(state);

            // Ok, now we have an inlet, check the values are right.
            REQUIRE(file->GetFilePath().filename() == "iolet.txt");
            REQUIRE(Approx(78.0) == converter.ConvertPressureToPhysicalUnits(file->GetPressureMin()));
            REQUIRE(Approx(82.0) == converter.ConvertPressureToPhysicalUnits(file->GetPressureMax()));
            auto expected = ApproxVector<PhysicalPosition>{-1.66017717834e-05, -4.58437586355e-05, -0.05};
            PhysicalPosition actual = converter.ConvertPositionToPhysicalUnits(file->GetPosition());
            REQUIRE(expected == actual);
            REQUIRE(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0) == file->GetNormal());

            // Set some target values for the density at various times.
            double targetStartDensity = converter.ConvertPressureToLatticeUnits(78.0) / Cs2;
            double targetMidDensity = converter.ConvertPressureToLatticeUnits(82.0) / Cs2;

            REQUIRE(Approx(targetStartDensity) == file->GetDensityMin());
            REQUIRE(Approx(targetStartDensity) == file->GetDensity(0));
            REQUIRE(Approx(targetMidDensity) == file->GetDensity(state.GetEndTimeStep() / 2));
        }

        SECTION("TestParabolicVelocityConstruct") {
            // Bootstrap ourselves an inoutlet, by loading config.xml.
            UncheckedSimConfigReader reader(Resource("config-velocity-iolet.xml").Path());
            auto config = reader.Read();
            configuration::SimBuilder builder(config);
            auto& converter = *builder.GetUnitConverter();

            auto p_vel_config = config.GetInlets()[0];
            auto p_vel = util::clone_dynamic_cast<lb::InOutLetParabolicVelocity>(builder.BuildIolet(p_vel_config));
            REQUIRE(p_vel != nullptr);

            // Bootstrap ourselves a unit converter, which the cosine needs in initialisation
            lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                            config.GetTotalTimeSteps());
            // at this stage, Initialise() has not been called, so the
            // unit converter will be invalid, so we will not be able to
            // convert to physical units.
            p_vel->Initialise(&converter);
            p_vel->Reset(state);

            // Check the IOLET contains the values expected given the file.
            /*
             * <inlet>
             <velocity radius="0.005" maximum="0.10"/>
             <normal x="0.0" y="0.0" z="1.0" />
             <position x="-1.66017717834e-05" y="-4.58437586355e-05" z="-0.05" />
             </inlet>
            */
            REQUIRE(5.0 == p_vel->GetRadius());
            REQUIRE(0.01 == p_vel->GetMaxSpeed());
            auto expected = ApproxVector<PhysicalPosition>{-1.66017717834e-05, -4.58437586355e-05, -0.05};
            REQUIRE(expected == converter.ConvertPositionToPhysicalUnits(p_vel->GetPosition()));
            REQUIRE(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0) == p_vel->GetNormal());
        }

        SECTION("TestWomersleyVelocityConstruct") {
            // Bootstrap ourselves a in inoutlet, by loading config.xml.
            UncheckedSimConfigReader reader(Resource("config_new_velocity_inlets.xml").Path());
            auto config = reader.Read();
            configuration::SimBuilder builder(config);
            auto& converter = *builder.GetUnitConverter();

            auto womersVelConfig = config.GetInlets()[0];
            auto womersVel = util::clone_dynamic_cast<lb::InOutLetWomersleyVelocity>(builder.BuildIolet(womersVelConfig));

            // Bootstrap ourselves a unit converter, which the cosine needs in initialisation
            lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                            config.GetTotalTimeSteps());
            double voxelSize = config.GetVoxelSize();
            // at this stage, Initialise() has not been called, so the
            // unit converter will be invalid, so we will not be able to
            // convert to physical units.
            womersVel->Initialise(&converter);
            womersVel->Reset(state);

            // Check the IOLET contains the values expected given the file.
            REQUIRE(10.0 == womersVel->GetRadius());
            REQUIRE(Approx(1e-6) == womersVel->GetPressureGradientAmplitude());
            REQUIRE(5.0 == womersVel->GetPeriod());
            REQUIRE(2.0 == womersVel->GetWomersleyNumber());
            REQUIRE(ApproxVector<PhysicalPosition>{0, 0, -0.05} ==
                    converter.ConvertPositionToPhysicalUnits(womersVel->GetPosition()));
            REQUIRE(ApproxVector<Dimensionless>{0.0, 0.0, 1.0} ==
                    womersVel->GetNormal());

            // Test that the analytical solution at r=R is 0
            PhysicalPosition pointAtCylinderWall(womersVel->GetRadius() * voxelSize, 0, -0.05);
            LatticePosition pointAtCylinderWallLatticeUnits(converter.ConvertPositionToLatticeUnits(pointAtCylinderWall));
            LatticeVelocity zeroVelAtWall(womersVel->GetVelocity(pointAtCylinderWallLatticeUnits,
                                                                 0));
            REQUIRE(ApproxVector<LatticeVelocity>{0, 0, 0}.Margin(1e-9) == zeroVelAtWall);

            // With a small enough Womersley number, the solution should
            // match the Poiseuille solution for the same pressure
            // difference after pi/2 radians and have changed direction
            // after 3*pi/2 radians.
            Dimensionless eta(1), nu(1);

            Dimensionless alpha = 1e-4;
            womersVel->SetWomersleyNumber(alpha);
            womersVel->SetPeriod(2 * PI * pow(womersVel->GetRadius(), 2) / (alpha * alpha * nu));

            LatticeSpeed poiseuilleSolution = womersVel->GetPressureGradientAmplitude()
                                              * pow(womersVel->GetRadius(), 2) / (4 * eta);

            LatticePosition pointAtCentrelineLatticeUnits(womersVel->GetPosition());

            {
                LatticeVelocity poiseuilleVelAtCentreLine(womersVel->GetVelocity(pointAtCentrelineLatticeUnits,
                                                                                 0.25 * womersVel->GetPeriod()));
                REQUIRE(ApproxVector<LatticeVelocity>{0.0, 0.0, poiseuilleSolution} == poiseuilleVelAtCentreLine);
            }

            {
                LatticeVelocity poiseuilleVelAtCentreLine(womersVel->GetVelocity(pointAtCentrelineLatticeUnits,
                                                                                 0.75 * womersVel->GetPeriod()));
                REQUIRE(ApproxVector<LatticeVelocity>{0, 0, -poiseuilleSolution} == poiseuilleVelAtCentreLine);
            }

        }

        SECTION("TestFileVelocityConstruct") {
            // We have to move to a tempdir, as the path specified in the
            // xml file is a relative path
            CopyResourceToTempdir("velocity_inlet.txt");
            MoveToTempdir();

            // Bootstrap ourselves a file velocity inlet, by loading an appropriate config file.
            UncheckedSimConfigReader reader(Resource("config_file_velocity_inlet.xml").Path());
            auto config = reader.Read();
            configuration::SimBuilder builder(config);
            auto& converter = *builder.GetUnitConverter();

            auto fileVelConfig = config.GetInlets()[0];
            auto fileVel = util::clone_dynamic_cast<lb::InOutLetFileVelocity>(builder.BuildIolet(fileVelConfig));

            lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                            config.GetTotalTimeSteps());
            // at this stage, Initialise() has not been called, so the
            // unit converter will be invalid, so we will not be able to
            // convert to physical units.
            fileVel->Initialise(&converter);
            fileVel->Reset(state);

            int FilePath_length = fileVel->GetFilePath().length();
            REQUIRE(fileVel->GetFilePath().substr(FilePath_length-19) == std::string("/velocity_inlet.txt"));
            REQUIRE(fileVel->GetRadius() == 20.0);
            auto expected = ApproxVector<PhysicalPosition>{0, 0, -0.05};
            PhysicalPosition actual = converter.ConvertPositionToPhysicalUnits(fileVel->GetPosition());
            REQUIRE(expected == actual);
            REQUIRE(ApproxVector<Dimensionless>{0.0, 0.0, 1.0} == fileVel->GetNormal());

            LatticePosition pointAtCentrelineLatticeUnits(fileVel->GetPosition());

            {
                LatticeVelocity velAtCentreLine(fileVel->GetVelocity(pointAtCentrelineLatticeUnits,
                                                                     converter.ConvertTimeToLatticeUnits(3.0)));
                PhysicalVelocity physVelCentreLine =
                        converter.ConvertVelocityToPhysicalUnits(velAtCentreLine);
                REQUIRE(ApproxVector<PhysicalVelocity>{0, 0, 0.01}.Margin(1e-9) == physVelCentreLine);
            }

            // Point equidistant to the centreline and wall
            LatticePosition pointEquidistant = fileVel->GetPosition();
            pointEquidistant[0] -= fileVel->GetRadius() / 2.0;

            {
                LatticeVelocity velAtPointEquidistant(fileVel->GetVelocity(pointEquidistant,
                                                                           converter.ConvertTimeToLatticeUnits(3.0)));
                PhysicalVelocity physVelPointEqui =
                        converter.ConvertVelocityToPhysicalUnits(velAtPointEquidistant);
                REQUIRE(ApproxVector<PhysicalVelocity>{0.0, 0.0, 0.0075}.Margin(1e-9) == physVelPointEqui);
            }

        }

        SECTION("TestIoletCoordinates") {
            // unit converter - make physical and lattice units the same
            util::UnitConverter units(1, 1, PhysicalPosition::Zero(), DEFAULT_FLUID_DENSITY_Kg_per_m3, 0.0);

            ConcreteIolet iolet;
            // normal
            util::Vector3D<Dimensionless> n(5, 7, -4);
            n.Normalise();
            iolet.SetNormal(n);
            // position
            PhysicalPosition c(7.77438796, 9.21293516, 9.87122463);
            iolet.SetPosition(units.ConvertPositionToLatticeUnits(c));
            iolet.Initialise(&units);
            lb::IoletExtraData extra(iolet);
            iolet.SetExtraData(&extra);

            // Convert the centre to iolet coords
            LatticePosition tmp = extra.WorldToIolet(c);
            // This should be zero
            REQUIRE(ApproxVector<LatticePosition>{}.Margin(1e-9) == tmp);

            // Make a point 3 lattice units along the normal.
            LatticePosition zEqThree = c + n * 3.0;
            tmp = extra.WorldToIolet(zEqThree);
            // This should be zero in x & y but 3 in z
            REQUIRE(ApproxVector<LatticePosition>{0,0,3}.Margin(1e-9) == tmp);
        }

    }

}

