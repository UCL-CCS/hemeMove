#ifndef HEMELB_REPORTING_REPORTER_H
#define HEMELB_REPORTING_REPORTER_H

#include <string>
#include <stdarg.h>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "Timers.h"
#include "Policies.h"
#include "lb/IncompressibilityChecker.h"
namespace hemelb
{
  namespace reporting
  {
    /**
     * Report generator class.
     * Class defining the creation of a report, intended for long-term archiving, describing what happened during a HemeLB run.
     * Accepts three policies as template arguments, defining:
     * @tparam ClockPolicy Clock to use in timing the performance
     * @tparam WriterPolicy How to write and format the report file.
     * @tparam CommsPolicy How to gather timing information across multiple processes
     * @tparam BroadcastPolicy The way we broadcast information across multiple processes
     * @todo CommsPolicy and BroadcastPolicy should be unified
     */
    template<class ClockPolicy, class WriterPolicy, class CommsPolicy, class BroadcastPolicy> class ReporterBase : public WriterPolicy,
                                                                                                                   public CommsPolicy
    {
      public:
        /**
         * Build a reporter capable of reporting on timings and other aspects of a HemeLB run.
         * @param path Path to write the report file to.
         * @param inputFile Input XML config file used to initialise the run.
         * @param aSiteCount Total count of sites used in the simulation.
         * @param timers Reference to list of timers used to measure performance.
         * @param aState Reference to state of ongoing simulation.
         */
        ReporterBase(const std::string &path,
                     const std::string &inputFile,
                     const long int aSiteCount,
                     const TimersBase<ClockPolicy, CommsPolicy>& timers,
                     const lb::SimulationState & aState,
                     const lb::IncompressibilityChecker<BroadcastPolicy>& aChecker);
        void Image(); //! Inform the reporter that an image has been saved.
        void Snapshot(); //! Inform the reporter that a simulation snapshot has been taken.
        void Write(); //! Write the report to disk, (or wherever the WriterPolicy decides.)
        void Stability(bool astability) //! Tell the reporter the current simulation stability state.
        {
          stability = astability;
        }
      private:
        bool doIo; //! Is this the processor which should write the report.
        unsigned int snapshotCount; //! Number of snapshots taken.
        unsigned int imageCount; //! Number of images written.
        long int siteCount; //! Total number of sites.
        bool stability; //! Stability of the simulation.
        const TimersBase<ClockPolicy, CommsPolicy> &timings; //! Reference to list of timers used to measure performance.
        const lb::SimulationState & state; //! Reference to state of ongoing simulation.
        const lb::IncompressibilityChecker<BroadcastPolicy>& incompressibilityChecker;
    };

    /**
     * Concrete realisation of the reporter with appropriate policies to be used.
     */
    typedef ReporterBase<HemeLBClockPolicy, FileWriterPolicy, MPICommsPolicy,
        net::PhasedBroadcastRegular<> > Reporter;
  }
}

#endif // HEMELB_REPORTING_REPORTER_H
