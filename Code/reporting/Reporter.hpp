#ifndef HEMELB_REPORTING_REPORTER_HPP
#define HEMELB_REPORTING_REPORTER_HPP
#include "Reporter.h"
namespace hemelb
{
  namespace reporting
  {
    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> ReporterBase<TimersPolicy,
        WriterPolicy, CommsPolicy>::ReporterBase(const std::string &name,
                                                 const std::string &inputFile,
                                                 const long int aSiteCount,
                                                 const TimersPolicy& timers,
                                                 const lb::SimulationState &aState) :
      WriterPolicy(name), snapshotCount(0), imageCount(0), siteCount(aSiteCount), stability(true),
          timings(timers), state(aState)
    {
      WriterPolicy::Print("***********************************************************\n");
      WriterPolicy::Print("Opening config file:\n %s\n", inputFile.c_str());
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<
        TimersPolicy, WriterPolicy, CommsPolicy>::Image()
    {
      imageCount++;
      WriterPolicy::Print("Image written: %u\n", imageCount);
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<
        TimersPolicy, WriterPolicy, CommsPolicy>::Snapshot()
    {
      snapshotCount++;
      WriterPolicy::Print("Snapshot written: %u\n", snapshotCount);
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<
        TimersPolicy, WriterPolicy, CommsPolicy>::Write()
    {

      // Note that CycleId is 1-indexed and will have just been incremented when we finish.
      unsigned long cycles = state.GetCycleId() - 1;

      WriterPolicy::Print("\n");
      WriterPolicy::Print("threads: %i, machines checked: %i\n\n",
                          CommsPolicy::GetProcessorCount(),
                          CommsPolicy::GetMachineCount());
      WriterPolicy::Print("topology depths checked: %i\n\n", CommsPolicy::GetDepths());
      WriterPolicy::Print("fluid sites: %li\n\n", siteCount);
      WriterPolicy::Print("cycles and total time steps: %lu, %lu \n\n",
                          cycles,
                          state.GetTimeStepsPassed() - 1);
      WriterPolicy::Print("time steps per second: %.3f\n\n",
                          (double) (state.GetTimeStepsPassed() - 1)
                              / timings[TimersPolicy::simulation].Get());

      if (!stability)
      {
        WriterPolicy::Print("Attention: simulation unstable with %lu timesteps/cycle\n",
                            state.GetTimeStepsPerCycle());
        WriterPolicy::Print("Simulation terminated\n");
      }

      WriterPolicy::Print("time steps per cycle: %lu\n", state.GetTimeStepsPerCycle());
      WriterPolicy::Print("\n");

      WriterPolicy::Print("\n");

      WriterPolicy::Print("\n");
      WriterPolicy::Print("total time (s):                            %.3f\n\n",
                           (timings[TimersPolicy::total].Get()));

      WriterPolicy::Print("Sub-domains info:\n\n");

      for (hemelb::proc_t n = 0; n < CommsPolicy::GetProcessorCount(); n++)
      {
        WriterPolicy::Print("rank: %lu, fluid sites: %lu\n",
                            (unsigned long) n,
                            (unsigned long) CommsPolicy::FluidSitesOnProcessor(n));
      }

      double normalisations[TimersPolicy::numberOfTimers] = { 1.0,
                                                              1.0,
                                                              1.0,
                                                              1.0,
                                                              (double) cycles,
                                                              (double) imageCount,
                                                              (double) cycles,
                                                              (double) cycles,
                                                              (double) snapshotCount,
                                                              1.0 };

      WriterPolicy::Print("\n\nPer-proc timing data (secs per [simulation,simulation,simulation,simulation,cycle,image,cycle,cycle,snapshot,simulation]): \n\n");
      WriterPolicy::Print("\t\tLocal \tMin \tMean \tMax\n");
      for (unsigned int ii = 0; ii < TimersPolicy::numberOfTimers; ii++)
      {
        WriterPolicy::Print("%s\t\t%.3g\t%.3g\t%.3g\t%.3g\n",
                            timerNames[ii].c_str(),
                            timings[ii].Get(),
                            timings.Mins()[ii] / normalisations[ii],
                            timings.Means()[ii] / normalisations[ii],
                            timings.Maxes()[ii] / normalisations[ii]);
      }
    }
  }
}
#endif //HEMELB_REPORTING_REPORTER_HPP
