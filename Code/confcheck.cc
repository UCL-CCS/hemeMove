// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <memory>
#include <string>

#include "configuration/SimConfig.h"
#include "io/ensure_hexfloat.h"
#include "net/MpiEnvironment.h"

using SimConfig = hemelb::configuration::SimConfig;

const char* usage = "confcheck hemelb_input_config.xml";

int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cerr << "Wrong number of arguments\n"
	      << "Usage: " << usage << std::endl;
    return 1;
  }

  const auto xml_path = std::string{argv[1]};

  hemelb::io::GlobalHexFloatLocale ensure_hexfloat;

  // When #755 is closed, remove MPI
  hemelb::net::MpiEnvironment mpi(argc, argv);
  try {
    auto conf = SimConfig::New(xml_path);
    return 0;
  } catch (std::exception& e) {
    std::cerr << "Error reading HemeLB configuration XML file: "
	      << e.what()
	      << std::endl;
    return 1;
  }
}
