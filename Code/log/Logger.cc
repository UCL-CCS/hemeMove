// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstdarg>
#include <sys/time.h>
#include <sys/resource.h>

#include "util/clock.h"
#include "net/mpi.h"
#include "log/Logger.h"

namespace hemelb::log
{
    // Use negative value to indicate uninitialised.
    int Logger::thisRank = -1;
    double Logger::startTime = -1.0;

    void Logger::Init()
    {
      // If Logger uninitialised
      if (thisRank < 0)
      {
        // Check that MPI is ready
        if (net::MpiEnvironment::Initialized())
        {
          thisRank = net::MpiCommunicator::World().Rank();
        }
        startTime = util::clock();
      }
    }

    template<>
    void Logger::LogInternal<OnePerCore>(std::string_view format, std::va_list args)
    {
      std::stringstream output;

      // Set the fill digit to be 0, so the integer 1 renders as 0000001
      output.fill('0');

      output << "[Rank " << std::setw(7) << thisRank << ", " << std::setiosflags(std::ios::fixed)
          << std::setprecision(1) << (util::clock() - startTime) << "s";

#ifdef HAVE_RUSAGE
      rusage usage;
      getrusage(RUSAGE_SELF, &usage);

      output << ", mem: " << std::setw(7) << usage.ru_maxrss;
#endif

      output << "]: " << format << '\n';

      std::string overFormat(output.str());

      std::vprintf(overFormat.c_str(), args);
    }

    template<>
    void Logger::LogInternal<Singleton>(std::string_view format, std::va_list args)
    {
        if (thisRank == 0)
        {
            std::printf("![%.1fs] ", util::clock() - startTime);
            std::vprintf(format.data(), args);
            std::printf("\n");
        }
    }
}
