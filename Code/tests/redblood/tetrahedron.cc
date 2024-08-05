// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <vector>
#include <sstream>
#include <string>
#include <iterator>
#include <exception>

#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "util/UnitConverter.h"

int get_depth(std::vector<std::string> const &_args)
{
  std::vector<std::string>::const_iterator i_found = std::find(_args.begin(),
                                                               _args.end(),
                                                               std::string("--depth"));

  if (i_found == _args.end())
  {
    i_found = std::find(_args.begin(), _args.end(), std::string("-d"));
  }

  if (i_found == _args.end())
  {
    return 0;
  }

  if ( (++i_found) == _args.end())
  {
    std::cerr << "Incorrect command-line argument";
    throw std::exception();
  }

  std::istringstream sstr(i_found->c_str());
  int result;
  sstr >> result;
  return result;
}

std::string get_output(std::vector<std::string> const &_args)
{
  std::vector<std::string>::const_iterator i_current = _args.begin();
  std::vector<std::string>::const_iterator const i_end = _args.end();

  for (; i_current != i_end; ++i_current)
    if ( (*i_current)[0] != '-')
    {
      return *i_current;
    }
    else if ( (++i_current) == i_end)
    {
      break;
    }

  return "";
}

bool pop_option(std::vector<std::string> &_args, std::string const &_option)
{
  std::vector<std::string>::const_iterator i_found = std::find(_args.begin(), _args.end(), _option);

  if (i_found == _args.end())
  {
    return false;
  }

  _args.erase(i_found);
  return true;
}

// Usage is $ tetrahedron [filename] [-d, --depth positive integer] [--vtk]
int main(int argc, char const *argv[])
{
  using namespace hemelb;
  using namespace hemelb::redblood;

  std::vector<std::string> args;
  std::copy(argv + 1, argv + argc, std::back_inserter(args));

  bool const dovtk = pop_option(args, "--vtk");
  int const depth(get_depth(args));
  std::string const output(get_output(args));

  Mesh result = tetrahedron(depth);

  auto io = [&]() -> std::unique_ptr<MeshIO> {
    if (dovtk) {
      return std::make_unique<VTKMeshIO>();
    }
    return std::make_unique<KruegerMeshIO>();
  }();

  if (output != "") {
    io->writeFile(output, *result.GetData(), nullptr);
  } else {
    std::cout << io->writeString(*result.GetData(), nullptr);
  }

  return 0;
}
