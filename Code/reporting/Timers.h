// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_TIMERS_H
#define HEMELB_REPORTING_TIMERS_H

#include <array>
#include <utility>

#include "reporting/Reportable.h"
#include "util/numerical.h"
#include "reporting/Policies.h"
#include "reporting/timers_fwd.h"

namespace hemelb::reporting
{

    template <typename T, typename Result>
    concept timer_function = std::is_default_constructible_v<T> && std::same_as<std::invoke_result_t<T>, Result>;

    /**
     * Timer which manages performance measurement for a single aspect of the code
     * @tparam ClockPolicy Policy defining how to get the current time
     */
    template<typename Timer>
    class TimerBase
    {
        // Assert rather than constrain type parameter to avoid
        // having concept in forward declaration.
        static_assert(timer_function<Timer, double>);

        Timer timer;
        double start = 0; //! Time when the timer was last started.
        double time = 0; //! Current running total time.

    public:
        /**
         * Starts with the timer stopped.
         */
        TimerBase() = default;

        /**
         * Get the current total time spent on this timer
         * @return current total time spent on this timer
         */
        double Get() const
        {
            return time;
        }
        /**
         * Force the current total time to an arbitrary value.
         * Use to reset if the simulation restarts
         * @param t The time to which to set the timer
         */
        void Set(double t)
        {
            time = t;
        }
        /**
         * Start the timer.
         */
        void Start()
        {
            start = timer();
        }
        /**
         * Stop the timer.
         */
        void Stop()
        {
            time += timer() - start;
        }
    };

    /**
     * Class which manages a set of timers timing aspects of a HemeLB run
     * @tparam ClockPolicy How to get the current time
     */
    template<class ClockPolicy>
    class TimersBase : public Reportable
    {
    public:
        class Timer : public TimerBase<ClockPolicy> {
            std::size_t index;
            std::string description;
            std::string name;
            friend TimersBase;
        public:
            auto Description() const {
                return description;
            }
        };
        /**
         * The set of possible timers
         */

        enum class NAMES
        {
            total = 0, //!< Total time
            initialDecomposition, //!< Initial seed decomposition
            domainDecomposition, //!< Time spent in parmetis domain decomposition
            fileRead, //!< Time spent in reading the geometry description file
            reRead, //!< Time spend in re-reading the geometry after second decomposition
            unzip, //!< Time spend in un-zipping
            moves, //!< Time spent moving things around post-parmetis
            parmetis, //!< Time spent in Parmetis
            latDatInitialise, //!< Time spent initialising the lattice data
            lb, //!< Time spent doing the core lattice boltzman simulation
            lb_calc, //!< Time spent doing calculations in the core lattice boltzmann simulation
            monitoring, //!< Time spent monitoring for stability, compressibility, etc.
            mpiSend, //!< Time spent sending MPI data
            mpiWait, //!< Time spent waiting for MPI
            simulation, //!< Total time for running the simulation,
            readNet,
            readParse,
            readBlock,
            readBlocksPrelim,
            readBlocksAll,
            moveForcingNumbers,
            moveForcingData,
            blockRequirements,
            moveCountsSending,
            moveDataSending,
            PopulateOptimisationMovesList,
            InitialGeometryRead,
            colloidInitialisation,
            colloidCommunicatePositions,
            colloidCommunicateVelocities,
            colloidCalculateForces,
            colloidUpdateCalculations,
            colloidOutput,
            extractionWriting,
            cellInitialisation,
            cellInsertion,
            computeNodeDistributions,
            exchangeCells,
            computeAndPostVelocities,
            receiveVelocitiesAndUpdate,
            updateDNC,
            computeAndPostForces,
            receiveForcesAndUpdate,
            updateCellAndWallInteractions,
            cellRemoval,
            cellListeners,
            graphComm,
            writeCheckpoint,
            last
            //!< last, this has to be the last element of the enumeration so it can be used to track cardinality
        };
        static constexpr std::size_t numberOfTimers = std::underlying_type_t<NAMES>(NAMES::last);
    private:
        int n_processes = 0;
        std::array<Timer, numberOfTimers> timers; //! The set of timers
        using DArray = std::array<double, numberOfTimers>;
        DArray maxes; //! Max across processes
        DArray mins; //! Min across processes
        DArray means; //! Average across processes

    public:
        TimersBase() {
#define HLB_TIMER_INIT(nm, descr) do { \
    auto i = std::size_t(NAMES::nm); \
    this->timers[i].index = i; \
    this->timers[i].name = #nm; \
    this->timers[i].description = descr; \
} while(0)

            HLB_TIMER_INIT(total, "Total");
            HLB_TIMER_INIT(initialDecomposition, "Seed Decomposition");
            HLB_TIMER_INIT(domainDecomposition, "Domain Decomposition");
            HLB_TIMER_INIT(fileRead, "File Read");
            HLB_TIMER_INIT(reRead, "Re Read");
            HLB_TIMER_INIT(unzip, "Unzip");
            HLB_TIMER_INIT(moves, "Moves");
            HLB_TIMER_INIT(parmetis, "Parmetis");
            HLB_TIMER_INIT(latDatInitialise, "Lattice Data initialisation");
            HLB_TIMER_INIT(lb, "Lattice Boltzmann");
            HLB_TIMER_INIT(lb_calc, "LB calc only");
            HLB_TIMER_INIT(monitoring, "Monitoring");
            HLB_TIMER_INIT(mpiSend, "MPI Send");
            HLB_TIMER_INIT(mpiWait, "MPI Wait");
            HLB_TIMER_INIT(simulation, "Simulation total");
            HLB_TIMER_INIT(readNet, "Reading communications");
            HLB_TIMER_INIT(readParse, "Parsing");
            HLB_TIMER_INIT(readBlock, "Read IO");
            HLB_TIMER_INIT(readBlocksPrelim, "Read Blocks prelim");
            HLB_TIMER_INIT(readBlocksAll, "Read blocks all");
            HLB_TIMER_INIT(moveForcingNumbers, "Move Forcing Counts");
            HLB_TIMER_INIT(moveForcingData, "Move Forcing Data");
            HLB_TIMER_INIT(blockRequirements, "Block Requirements");
            HLB_TIMER_INIT(moveCountsSending, "Move Counts Sending");
            HLB_TIMER_INIT(moveDataSending, "Move Data Sending");
            HLB_TIMER_INIT(PopulateOptimisationMovesList, "Populating moves list for decomposition optimisation");
            HLB_TIMER_INIT(InitialGeometryRead, "Initial geometry reading");
            HLB_TIMER_INIT(colloidInitialisation, "Colloid initialisation");
            HLB_TIMER_INIT(colloidCommunicatePositions, "Colloid position communication");
            HLB_TIMER_INIT(colloidCommunicateVelocities, "Colloid velocity communication");
            HLB_TIMER_INIT(colloidCalculateForces, "Colloid force calculations");
            HLB_TIMER_INIT(colloidUpdateCalculations, "Colloid calculations for updating");
            HLB_TIMER_INIT(colloidOutput, "Colloid outputting");
            HLB_TIMER_INIT(extractionWriting, "Extraction writing");
            HLB_TIMER_INIT(cellInitialisation, "RBC initialisation");
            HLB_TIMER_INIT(cellInsertion, "RBC insertion");
            HLB_TIMER_INIT(computeNodeDistributions, "Compute node distributions");
            HLB_TIMER_INIT(exchangeCells, "Exchange cells");
            HLB_TIMER_INIT(computeAndPostVelocities, "Compute local velocities and post them");
            HLB_TIMER_INIT(receiveVelocitiesAndUpdate, "Receive velocities and update non local contributions");
            HLB_TIMER_INIT(updateDNC, "Update divide and conquer");
            HLB_TIMER_INIT(computeAndPostForces, "Compute local forces and post them");
            HLB_TIMER_INIT(receiveForcesAndUpdate, "Receive forces and update non local contributions");
            HLB_TIMER_INIT(updateCellAndWallInteractions, "Update cell-cell and cell-wall interactions");
            HLB_TIMER_INIT(cellRemoval, "Remove cells");
            HLB_TIMER_INIT(cellListeners, "Notify cell listeners");
            HLB_TIMER_INIT(graphComm, "Create graph communicator");
            HLB_TIMER_INIT(writeCheckpoint, "Write checkpoint(s)");
#undef HLB_TIMER_INIT

            for (auto const& t: timers) {
                if (t.description.empty())
                    throw (Exception() << "Missing description for timer");
            }
        }
        ~TimersBase() noexcept override = default;

#define HLB_TIMER_ACCESSOR(name) Timer& name() { \
        return timers[std::size_t(NAMES::name)];       \
        }
        HLB_TIMER_ACCESSOR(total)
        HLB_TIMER_ACCESSOR(initialDecomposition)
        HLB_TIMER_ACCESSOR(domainDecomposition)
        HLB_TIMER_ACCESSOR(fileRead)
        HLB_TIMER_ACCESSOR(reRead)
        HLB_TIMER_ACCESSOR(unzip)
        HLB_TIMER_ACCESSOR(moves)
        HLB_TIMER_ACCESSOR(parmetis)
        HLB_TIMER_ACCESSOR(latDatInitialise)
        HLB_TIMER_ACCESSOR(lb)
        HLB_TIMER_ACCESSOR(lb_calc)
        HLB_TIMER_ACCESSOR(monitoring)
        HLB_TIMER_ACCESSOR(mpiSend)
        HLB_TIMER_ACCESSOR(mpiWait)
        HLB_TIMER_ACCESSOR(simulation)
        HLB_TIMER_ACCESSOR(readNet)
        HLB_TIMER_ACCESSOR(readParse)
        HLB_TIMER_ACCESSOR(readBlock)
        HLB_TIMER_ACCESSOR(readBlocksPrelim)
        HLB_TIMER_ACCESSOR(readBlocksAll)
        HLB_TIMER_ACCESSOR(moveForcingNumbers)
        HLB_TIMER_ACCESSOR(moveForcingData)
        HLB_TIMER_ACCESSOR(blockRequirements)
        HLB_TIMER_ACCESSOR(moveCountsSending)
        HLB_TIMER_ACCESSOR(moveDataSending)
        HLB_TIMER_ACCESSOR(PopulateOptimisationMovesList)
        HLB_TIMER_ACCESSOR(InitialGeometryRead)
        HLB_TIMER_ACCESSOR(colloidInitialisation)
        HLB_TIMER_ACCESSOR(colloidCommunicatePositions)
        HLB_TIMER_ACCESSOR(colloidCommunicateVelocities)
        HLB_TIMER_ACCESSOR(colloidCalculateForces)
        HLB_TIMER_ACCESSOR(colloidUpdateCalculations)
        HLB_TIMER_ACCESSOR(colloidOutput)
        HLB_TIMER_ACCESSOR(extractionWriting)
        HLB_TIMER_ACCESSOR(cellInitialisation)
        HLB_TIMER_ACCESSOR(cellInsertion)
        HLB_TIMER_ACCESSOR(computeNodeDistributions)
        HLB_TIMER_ACCESSOR(exchangeCells)
        HLB_TIMER_ACCESSOR(computeAndPostVelocities)
        HLB_TIMER_ACCESSOR(receiveVelocitiesAndUpdate)
        HLB_TIMER_ACCESSOR(updateDNC)
        HLB_TIMER_ACCESSOR(computeAndPostForces)
        HLB_TIMER_ACCESSOR(receiveForcesAndUpdate)
        HLB_TIMER_ACCESSOR(updateCellAndWallInteractions)
        HLB_TIMER_ACCESSOR(cellRemoval)
        HLB_TIMER_ACCESSOR(cellListeners)
        HLB_TIMER_ACCESSOR(graphComm)
        HLB_TIMER_ACCESSOR(writeCheckpoint)
#undef HLB_TIMER_ACCESSOR
        /**
         * Max across all processes.
         * Following the sharing of timing data between processes, the max time across all processes for each timer.
         * @return the max time across all processes for each timer.
         */
        auto& Maxes() const
        {
            return maxes;
        }
        /**
         * Min across all processes.
         * Following the sharing of timing data between processes, the minimum time across all processes for each timer.
         * @return the minimum time across all processes for each timer.
         */
        auto& Mins() const
        {
            return mins;
        }
        /**
         * Averages across all processes.
         * Following the sharing of timing data between processes, the average time across all processes for each timer.
         * @return the average time across all processes for each timer.
         */
        auto& Means() const
        {
            return means;
        }

        auto begin() {
            return timers.begin();
        }
        auto begin() const {
            return timers.begin();
        }
        auto end() {
            return timers.end();
        }
        auto end() const {
            return timers.end();
        }

        /**
         * Share timing information across timers
         */
        template <typename Communicator>
        void Reduce(Communicator&& comm);

        void Report(Dict& dictionary) override;

    private:
    };

    using Timer = TimerBase<HemeLBClockPolicy>;
    using Timers = TimersBase<HemeLBClockPolicy>;

}

#endif //HEMELB_REPORTING_TIMERS_H
