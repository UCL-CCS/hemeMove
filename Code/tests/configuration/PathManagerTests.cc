// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>

#include "configuration/CommandLine.h"
#include "configuration/PathManager.h"

#include "tests/helpers/FolderTestFixture.h"

namespace fs = std::filesystem;

namespace hemelb::tests
{

    std::unique_ptr<configuration::PathManager> setup(const char* confXml) {
        const int processorCount = 5;
        const int argc = 3;
        const char* argv[] = {
                "hemelb",
                "-in",
                confXml,
        };

        auto cl = configuration::CommandLine(argc, argv);
        return std::make_unique<configuration::PathManager>(cl, true, processorCount);
    }

    TEST_CASE_METHOD(helpers::FolderTestFixture, "PathManager") {
        SECTION("Local config XML path") {
            auto fileManager = setup("config.xml");

            REQUIRE(fs::exists("results"));

            REQUIRE(fs::absolute("results") == fileManager->GetReportPath());

        }

        SECTION("Abs path to config XML") {
            auto const targetConfig = GetTempdir() / "config.xml"; // note this resource doesn't exist -- not a problem
            ReturnToOrigin(); // even if we're not in current dir, explicit path should cause the results to be created in the tmpdir
            auto fileManager = setup(targetConfig.c_str());
            MoveToTempdir(); // go back to the tempdir and check the files were created in the right place

            REQUIRE(fs::exists("results"));
            REQUIRE(GetTempdir() / "results" == fileManager->GetReportPath());
        }

    }

}
