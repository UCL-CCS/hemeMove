// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LOG_LOGGER_H
#define HEMELB_LOG_LOGGER_H

#include <vector>
#include <cstdarg>
#include <sstream>
#include <string>

#include "build_info.h"

namespace hemelb::log
{
    enum LogLevel
    {
      //! Unrecoverable error
      Critical = 0,
      //! Recoverable error or unexpected state
      Error = 1,
      //! Undesirable but not necessarily wrong states
      Warning = 2,
      //! A conservative set of interesting run time events
      Info = 3,
      //! Detailed information on runtime flow
      //! This level should also be used for any non-essential computation.
      Debug = 4,
      //! Low-level trace information, including actual data.
      Trace = 5
    };

    namespace detail {
        consteval LogLevel default_log_level() {
            auto l = build_info::LOG_LEVEL;
            if (l == "Critical")
                return Critical;
            else if (l == "Error")
                return Error;
            else if (l == "Warning")
                return Warning;
            else if (l == "Info")
                return Info;
            else if (l == "Debug")
                return Debug;
            else if (l == "Trace")
                return Trace;
            else
                throw "Invalid choice for LOG_LEVEL";
        }
    }

    enum LogType
    {
      Singleton,
      OnePerCore
    };

    class Logger
    {
    public:
        template<LogLevel queryLogLevel>
        static constexpr bool ShouldDisplay()
        {
            return queryLogLevel <= currentLogLevel;
        }

        static void Init();

        template<LogLevel queryLogLevel, LogType logType>
        static void Log(std::string_view format, ...)
        {
          if constexpr (queryLogLevel <= currentLogLevel)
          {
            va_list args;
            va_start(args, format);
            LogInternal<logType>(format, args);
            va_end(args);
          }
        }

    private:
        template<LogType>
        static void LogInternal(std::string_view format, va_list args);

        static constexpr LogLevel currentLogLevel = detail::default_log_level();
        static int thisRank;
        static double startTime;
    };


    template<class... ARGS>
    void capture(std::string_view message, ARGS &&... args)
    {
        std::ostringstream sstr;
        (sstr << ... << args);
        Logger::Log<Debug, OnePerCore>(
                message,
                sstr.str().c_str()
        );
    }
}

//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE(VARIABLE) ::hemelb::log::capture(#VARIABLE " := %s", VARIABLE)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE2(V0, V1) ::hemelb::log::capture(#V0 " := %s " #V1 " := %s", V0, V1)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE3(V0, V1, V2)     \
  ::hemelb::log::capture(#V0 " := %s " #V1 " := %s " #V2 " := %s", V0, V1, V2)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE4(V0, V1, V2, V3) \
  ::hemelb::log::capture(#V0 " := %s " #V1 " := %s " #V2 " := %s " #V3 ": = %s", V0, V1, V2, V3)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE5(V0, V1, V2, V3, V4) \
  ::hemelb::log::capture(                     \
    #V0 " := %s " #V1 " := %s " #V2 " := %s " #V3 ": = %s " #V4 ": = %s", V0, V1, V2, V3, V4)

#endif /* HEMELB_LOG_LOGGER_H */
