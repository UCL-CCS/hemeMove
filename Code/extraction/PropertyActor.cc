// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/PropertyActor.h"
#include "reporting/Timers.h"

namespace hemelb::extraction
{
    PropertyActor::PropertyActor(std::shared_ptr<lb::SimulationState const> simState,
                                 const std::vector<PropertyOutputFile>& propertyOutputs,
                                 std::shared_ptr<IterableDataSource> dataSource, reporting::Timers& timers,
                                 net::IOCommunicator const& ioComms) :
        simulationState(std::move(simState)), timers(timers)
    {
      propertyWriter = std::make_unique<PropertyWriter>(std::move(dataSource), propertyOutputs, ioComms);
    }

    PropertyActor::~PropertyActor() = default;

    void PropertyActor::SetRequiredProperties(lb::MacroscopicPropertyCache& propertyCache)
    {
        auto& propertyOutputs = propertyWriter->GetPropertyOutputs();

        // Iterate over each property output spec.
        for (auto& propertyOutput : propertyOutputs)
        {
            // Only consider the ones that are being written this iteration.
            if (propertyOutput.ShouldWrite(simulationState->GetTimeStep()))
            {
                auto& outputFile = propertyOutput.GetOutputSpec();

                // Iterate over each field.
                for (auto&& fieldSpec: outputFile.fields)
                {
                    // Set the cache to calculate each required field.
                    overload_visit(
                            fieldSpec.src,
                            [&](source::Pressure) {
                                propertyCache.densityCache.SetRefreshFlag();
                            },
                            [&](source::Velocity) {
                                propertyCache.velocityCache.SetRefreshFlag();
                            },
                            [&](source::ShearStress) {
                                propertyCache.wallShearStressMagnitudeCache.SetRefreshFlag();
                            },
                            [&](source::VonMisesStress) {
                                propertyCache.vonMisesStressCache.SetRefreshFlag();
                            },
                            [&](source::ShearRate) {
                                propertyCache.shearRateCache.SetRefreshFlag();
                            },
                            [&](source::StressTensor) {
                                propertyCache.stressTensorCache.SetRefreshFlag();
                            },
                            [&](source::Traction) {
                                propertyCache.tractionCache.SetRefreshFlag();
                            },
                            [&](source::TangentialProjectionTraction) {
                                propertyCache.tangentialProjectionTractionCache.SetRefreshFlag();
                            },
                            [](source::Distributions) {
                                // We don't actually have to cache anything to get the distribution.
                            },
                            [](source::MpiRank) {
                                // We don't actually have to cache anything to get the rank.
                            }
                    );

                }
            }
        }
    }

    void PropertyActor::EndIteration()
    {
      timers.extractionWriting().Start();
      propertyWriter->Write(simulationState->GetTimeStep(), simulationState->GetEndTimeStep());
      timers.extractionWriting().Stop();
    }

}
