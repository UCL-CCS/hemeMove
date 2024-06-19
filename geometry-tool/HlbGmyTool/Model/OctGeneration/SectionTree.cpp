// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "SectionTree.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include "H5.h"
#include "enumerate.hpp"

namespace hemelb {
namespace H5 {

using namespace gmytool::oct;

template <>
DataTypeSharedPtr DataTypeTraits<SVector>::GetType() {
  return DataType::Array<float>({3});
}

template <>
DataTypeSharedPtr DataTypeTraits<std::array<Link, 26>>::GetType() {
  return DataType::Array<Link>({26});
}

template <>
DataTypeSharedPtr DataTypeTraits<Link>::GetType() {
  hid_t link_id = H5Tcreate(H5T_COMPOUND, sizeof(Link));
  H5Tinsert(link_id, "type", HOFFSET(Link, type), H5T_NATIVE_INT);
  H5Tinsert(link_id, "dist", HOFFSET(Link, dist), H5T_NATIVE_FLOAT);
  H5Tinsert(link_id, "id", HOFFSET(Link, id), H5T_NATIVE_INT);

  return DataTypeSharedPtr(new DataType(link_id));
}
}  // namespace H5

namespace gmytool::oct {

SectionTree::SectionTree(size_t nl)
    : nLevels(nl), indices(nl + 1), counts(nl + 1) {}

const SectionTree::Tree& SectionTree::GetTree() const {
  return indices;
}

auto SectionTree::FindIndex(Int i, Int j, Int k) const -> IndT {
  Int cur_level = nLevels;
  // Root level (==nLevel) has implicit offset of zero
  IndT cur_offset = 0;
  while (cur_level) {
    // Get the local index
    auto lIndex = LocalOffset(i, j, k, cur_level - 1);
    cur_offset = indices[cur_level][cur_offset + lIndex];
    if (cur_offset == NA()) {
      break;
    }
    --cur_level;
  }
  return cur_offset;
}

// H5 does not allow the chunk size to be bigger than the data set's
// size. Make a plist that satisfies this constraint.
H5::PListSharedPtr compression_plist(int size) {
  const int CHUNK = 512;
  if (size < CHUNK)
    // Don't bother with compression
    return H5::PList::Default();

  auto pl = H5::PList::DatasetCreate();
  pl->SetChunk({CHUNK});
  pl->SetDeflate(6);
  return pl;
}

template <class T>
void Section<T>::write(H5::GroupPtr grp) const {
  grp->CreateWriteDataSet("offsets", offsets,
                          compression_plist(offsets.size()));
  grp->CreateWriteDataSet("counts", counts, compression_plist(counts.size()));
  grp->CreateWriteDataSet("data", data, compression_plist(data.size()));
}

void SectionTree::Write(const std::string& fn) const {
  auto outfile = H5::File::Create(fn, H5F_ACC_TRUNC);
  auto root = outfile->CreateGroup("hemelb")->CreateGroup("geometry");
  root->SetAttribute("version", 1);
  root->SetAttribute("levels", int(nLevels));

  for (auto&& [i, level] : const_enumerate(indices)) {
    // Each level is a single dataset
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(4) << i;
    // Enable compression with a chunk size of 1k
    root->CreateWriteDataSet(ss.str(), level, compression_plist(level.size()));
  }

  wall_normals.write(root->CreateGroup("wall_normals"));
  links.write(root->CreateGroup("links"));
}

}  // namespace gmytool::oct
}  // namespace hemelb
