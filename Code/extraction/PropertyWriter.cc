// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/PropertyWriter.h"

namespace hemelb::extraction
{
    PropertyWriter::PropertyWriter(std::shared_ptr<IterableDataSource> dataSource,
                                   const std::vector<PropertyOutputFile>& propertyOutputs,
                                   net::IOCommunicator const& ioComms)
    {
      for (const auto& propertyOutput : propertyOutputs)
      {
        localPropertyOutputs.emplace_back(dataSource, propertyOutput, ioComms);
      }
    }

    const std::vector<LocalPropertyOutput>& PropertyWriter::GetPropertyOutputs() const
    {
      return localPropertyOutputs;
    }

    void PropertyWriter::Write(unsigned long iterationNumber, unsigned long totalSteps)
    {
      for (auto& localPropertyOutput : localPropertyOutputs)
      {
        localPropertyOutput.Write((uint64_t) iterationNumber, totalSteps);
      }
    }
}
