// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPIERROR_H
#define HEMELB_NET_MPIERROR_H

#include <functional>
#include <source_location>

#include <mpi.h>

#include "Exception.h"

namespace hemelb::net
{
    /**
     * Indicate an error to do with MPI.
     *
     * Will be thrown by the HEMELB_MPI_CALL macro, as in:
     *   HEMELB_MPI_CALL(MPI_Send, (buffer, count, dtype, toRank, tag, communicator) );
     *
     */
    class MpiError : public ::hemelb::Exception
    {
      public:
        MpiError(const char* mpiFunc_, int errorCode_, const char* fileName_, const int lineNo);

      private:
        const char* mpiFunc;
        const int errorCode;
        const char* fileName;
        const int lineNo;
    };

    /**
     * This helper will call the MPI_* function and check the return
     * code, ensuring that errors are handled.
     *
     * Use:
     * MpiCall{MPI_Bcast}(buffer, count, dtype, root, comm):
     */
    template <typename FuncT>
    struct MpiCall {
        FuncT* func;
        std::source_location loc;

        MpiCall(FuncT* f, std::source_location l = std::source_location::current())
            : func(std::move(f)), loc(l) {
        }

        template <typename... ArgTs>
        requires std::invocable<FuncT, ArgTs...>
        void operator()(ArgTs&&... args) {
            int res = std::invoke(func, std::forward<ArgTs>(args)...);
            if (res != MPI_SUCCESS)
              throw MpiError(loc.function_name(), res, loc.file_name(), loc.line());
        }
    };
}

// Macro version of the above - deprecated
#define HEMELB_MPI_CALL( mpiFunc, args ) do \
{ \
  int _check_result = mpiFunc args; \
  if (_check_result != MPI_SUCCESS) \
    throw ::hemelb::net::MpiError(#mpiFunc, _check_result, __FILE__, __LINE__); \
} while(0)

#endif // HEMELB_NET_MPIERROR_H
