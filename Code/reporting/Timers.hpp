// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_TIMERS_HPP
#define HEMELB_REPORTING_TIMERS_HPP

#include "reporting/Timers.h"

namespace hemelb::reporting
{
    template<class ClockPolicy>
    template <typename Communicator>
    void TimersBase<ClockPolicy>::Reduce(Communicator &&comm)
    {
        n_processes = comm.Size();
        DArray timings;
        for (unsigned int i = 0; i < numberOfTimers; i++)
            timings[i] = timers[i].Get();

        using span = std::span<double, numberOfTimers>;
        using cspan = std::span<double const, numberOfTimers>;
        comm.Reduce(span(maxes), cspan(timings), MPI_MAX, 0);
        comm.Reduce(span(means), cspan(timings), MPI_SUM, 0);
        comm.Reduce(span(mins), cspan(timings), MPI_MIN, 0);

        for (auto& m: means)
            m /= n_processes;
    }

    template<class ClockPolicy>
    void TimersBase<ClockPolicy>::Report(Dict& dictionary)
    {
      dictionary.SetIntValue("THREADS", n_processes);

      for (unsigned int ii = 0; ii < numberOfTimers; ii++)
      {
        Dict timer = dictionary.AddSectionDictionary("TIMER");
        timer.SetValue("NAME", timers[ii].description);
        timer.SetFormattedValue("LOCAL", "%.3g", timers[ii].Get());
        timer.SetFormattedValue("MIN", "%.3g", Mins()[ii]);
        timer.SetFormattedValue("MEAN", "%.3g", Means()[ii]);
        timer.SetFormattedValue("MAX", "%.3g", Maxes()[ii]);
      }
    }

}
#endif // HEMELB_REPORTING_TIMERS_HPP
