// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "H5.h"

#include <iostream>

#define H5_CLOSETHIS(func, nothr)            \
  if (m_Id != H5I_INVALID_HID) {             \
    herr_t ret = func(m_Id);                 \
    if (ret < 0)                             \
      if constexpr (nothr) {                 \
        std::cerr << H5_ERRORMSG(func) "\n"; \
        std::terminate();                    \
      } else {                               \
        throw Error(H5_ERRORMSG(func));      \
      }                                      \
    m_Id = H5I_INVALID_HID;                  \
  }

namespace hemelb::H5 {

Object::Object() : m_Id(H5I_INVALID_HID) {}
Object::Object(hid_t id) : m_Id(id) {}
Object::~Object() noexcept {}

PList::PList() : Object(H5P_DEFAULT) {}
PList::PList(hid_t cls) : Object() {
  H5_CONSTRUCT(m_Id, H5Pcreate, (cls));
}
PList::~PList() {
  H5_CLOSETHIS(H5Pclose, true);
}
void PList::Close() {
  H5_CLOSETHIS(H5Pclose, false);
}

/// Default options
PListSharedPtr PList::Default() {
  return std::make_shared<PList>();
}

/// Properties for object creation
PListSharedPtr PList::ObjectCreate() {
  return std::make_shared<PList>(H5P_OBJECT_CREATE);
}

/// Properties for file creation
PListSharedPtr PList::FileCreate() {
  return std::make_shared<PList>(H5P_FILE_CREATE);
}

/// Properties for file access
PListSharedPtr PList::FileAccess() {
  return std::make_shared<PList>(H5P_FILE_ACCESS);
}

/// Properties for dataset creation
PListSharedPtr PList::DatasetCreate() {
  return std::make_shared<PList>(H5P_DATASET_CREATE);
}

/// Properties for dataset access
PListSharedPtr PList::DatasetAccess() {
  return std::make_shared<PList>(H5P_DATASET_ACCESS);
}

/// Properties for raw data transfer
PListSharedPtr PList::DatasetXfer() {
  return std::make_shared<PList>(H5P_DATASET_XFER);
}

/// Properties for file mounting
PListSharedPtr PList::FileMount() {
  return std::make_shared<PList>(H5P_FILE_MOUNT);
}

/// Properties for group creation
PListSharedPtr PList::GroupCreate() {
  return std::make_shared<PList>(H5P_GROUP_CREATE);
}

/// Properties for group access
PListSharedPtr PList::GroupAccess() {
  return std::make_shared<PList>(H5P_GROUP_ACCESS);
}

/// Properties for datatype creation
PListSharedPtr PList::DatatypeCreate() {
  return std::make_shared<PList>(H5P_DATATYPE_CREATE);
}

/// Properties for datatype access
PListSharedPtr PList::DatatypeAccess() {
  return std::make_shared<PList>(H5P_DATATYPE_ACCESS);
}

/// Properties for character encoding when encoding strings or object names
PListSharedPtr PList::StringCreate() {
  return std::make_shared<PList>(H5P_STRING_CREATE);
}

/// Properties for attribute creation
PListSharedPtr PList::AttributeCreate() {
  return std::make_shared<PList>(H5P_ATTRIBUTE_CREATE);
}

/// Properties governing the object copying process
PListSharedPtr PList::ObjectCopy() {
  return std::make_shared<PList>(H5P_OBJECT_COPY);
}

/// Properties governing link creation
PListSharedPtr PList::LinkCreate() {
  return std::make_shared<PList>(H5P_LINK_CREATE);
}

/// Properties governing link traversal when accessing objects
PListSharedPtr PList::LinkAccess() {
  return std::make_shared<PList>(H5P_LINK_ACCESS);
}
void PList::SetChunk(const std::vector<hsize_t>& dims) {
  H5_CALL(H5Pset_chunk, (m_Id, dims.size(), &dims[0]));
}
void PList::SetDeflate(const unsigned level) {
  H5_CALL(H5Pset_deflate, (m_Id, level));
}

GroupSharedPtr CanHaveGroupsDataSets::CreateGroup(const std::string& name,
                                                  PListSharedPtr createPL,
                                                  PListSharedPtr accessPL) {
  hid_t grp;
  H5_CONSTRUCT(
      grp, H5Gcreate,
      (m_Id, name.c_str(), H5P_DEFAULT, createPL->GetId(), accessPL->GetId()));
  return std::make_shared<Group>(grp);
}

DataSetSharedPtr CanHaveGroupsDataSets::CreateDataSet(const std::string& name,
                                                      DataTypeSharedPtr type,
                                                      DataSpaceSharedPtr space,
                                                      PListSharedPtr createPL,
                                                      PListSharedPtr accessPL) {
  hid_t ds;
  H5_CONSTRUCT(ds, H5Dcreate,
               (m_Id, name.c_str(), type->GetId(), space->GetId(), H5P_DEFAULT,
                createPL->GetId(), accessPL->GetId()));

  return std::make_shared<DataSet>(ds);
}

// Open an existing group.
// The accessPL can be omitted to use the defaults
GroupSharedPtr CanHaveGroupsDataSets::OpenGroup(const std::string& name,
                                                PListSharedPtr accessPL) const {
  hid_t grp;
  H5_CONSTRUCT(grp, H5Gopen2, (m_Id, name.c_str(), accessPL->GetId()));
  return std::make_shared<Group>(grp);
}

// Open an existing dataset
// The accessPL can be omitted to use the defaults
DataSetSharedPtr CanHaveGroupsDataSets::OpenDataSet(
    const std::string& name,
    PListSharedPtr accessPL) const {
  hid_t ds;
  H5_CONSTRUCT(ds, H5Dopen2, (m_Id, name.c_str(), accessPL->GetId()));
  return std::make_shared<DataSet>(ds);
}

LinkIterator CanHaveGroupsDataSets::begin() {
  return {this, 0};
}

LinkIterator CanHaveGroupsDataSets::end() {
  return {this, GetNumElements()};
}

AttributeSharedPtr CanHaveAttributes::CreateAttribute(
    const std::string& name,
    DataTypeSharedPtr type,
    DataSpaceSharedPtr space) {
  return Attribute::Create(m_Id, name, type, space);
}

AttributeSharedPtr CanHaveAttributes::OpenAttribute(const std::string& name) {
  return Attribute::Open(m_Id, name);
}

hsize_t CanHaveAttributes::GetNumAttr() const {
  H5O_info_t info;
  H5_CALL(H5Oget_info, (m_Id, &info));
  return info.num_attrs;
}

AttrIterator CanHaveAttributes::attr_begin() const {
  return {this, 0};
}
AttrIterator CanHaveAttributes::attr_end() const {
  return {this, GetNumAttr()};
}

DataSpaceSharedPtr DataSpace::Null() {
  auto ans = std::make_shared<DataSpace>();
  H5_CONSTRUCT(ans->m_Id, H5Screate, (H5S_NULL));
  return ans;
}

DataSpaceSharedPtr DataSpace::Scalar() {
  auto ans = std::make_shared<DataSpace>();
  H5_CONSTRUCT(ans->m_Id, H5Screate, (H5S_SCALAR));
  return ans;
}

DataSpaceSharedPtr DataSpace::OneD(hsize_t size) {
  auto ans = std::make_shared<DataSpace>();
  H5_CONSTRUCT(ans->m_Id, H5Screate_simple, (1, &size, NULL));
  return ans;
}

DataSpace::DataSpace() : Object() {}

DataSpace::DataSpace(hid_t id) : Object(id) {}

DataSpace::DataSpace(const std::vector<hsize_t>& dims) : Object() {
  int rank = dims.size();
  H5_CONSTRUCT(m_Id, H5Screate_simple, (rank, &dims[0], NULL));
}

DataSpace::DataSpace(const hsize_t size, const hsize_t max) : Object() {
  const hsize_t* max_p = NULL;
  if (max != (H5S_UNLIMITED - 1))
    max_p = &max;
  H5_CONSTRUCT(m_Id, H5Screate_simple, (1, &size, max_p));
}

DataSpace::DataSpace(const std::vector<hsize_t>& dims,
                     const std::vector<hsize_t>& max_dims)
    : Object() {
  int rank = dims.size();
  H5_CONSTRUCT(m_Id, H5Screate_simple, (rank, &dims[0], &max_dims[0]));
}

DataSpace::~DataSpace() {
  H5_CLOSETHIS(H5Sclose, true);
}

void DataSpace::Close() {
  H5_CLOSETHIS(H5Sclose, false);
}

void DataSpace::SelectRange(const hsize_t start, const hsize_t count) {
  H5_CALL(H5Sselect_hyperslab,
          (m_Id, H5S_SELECT_SET, &start, NULL, &count, NULL));
}

DataType::DataType(hid_t id) : Object(id) {}
DataTypeSharedPtr DataType::String(size_t len) {
  DataTypeSharedPtr s1 = PredefinedDataType::CS1();
  DataTypeSharedPtr ans = s1->Copy();
  if (len == 0)
    len = H5T_VARIABLE;
  H5_CALL(H5Tset_size, (ans->GetId(), len));
  return ans;
}

DataType::~DataType() {
  H5_CLOSETHIS(H5Tclose, true);
}

void DataType::Close() {
  H5_CLOSETHIS(H5Tclose, false);
}

DataTypeSharedPtr DataType::Copy() const {
  hid_t ans_id = H5I_INVALID_HID;
  H5_CONSTRUCT(ans_id, H5Tcopy, (m_Id));
  return std::make_shared<DataType>(ans_id);
}

DataTypeSharedPtr PredefinedDataType::CS1() {
  return std::make_shared<PredefinedDataType>(H5T_C_S1);
}

PredefinedDataType::PredefinedDataType(hid_t id) : DataType(id) {}

void PredefinedDataType::Close() {
  // No-op
  m_Id = H5I_INVALID_HID;
}

PredefinedDataType::~PredefinedDataType() {
  // To ensure base class does not try to call H5Tclose
  m_Id = H5I_INVALID_HID;
}

template <>
const hid_t DataTypeTraits<char>::NativeType = H5T_NATIVE_CHAR;

template <>
const hid_t DataTypeTraits<int>::NativeType = H5T_NATIVE_INT;

template <>
const hid_t DataTypeTraits<unsigned int>::NativeType = H5T_NATIVE_UINT;

template <>
const hid_t DataTypeTraits<unsigned long>::NativeType = H5T_NATIVE_ULONG;

template <>
const hid_t DataTypeTraits<unsigned long long>::NativeType = H5T_NATIVE_ULLONG;

template <>
const hid_t DataTypeTraits<float>::NativeType = H5T_NATIVE_FLOAT;
template <>
const hid_t DataTypeTraits<double>::NativeType = H5T_NATIVE_DOUBLE;

Attribute::Attribute(hid_t id) : Object{id} {}

AttributeSharedPtr Attribute::Create(hid_t parent,
                                     const std::string& name,
                                     DataTypeSharedPtr type,
                                     DataSpaceSharedPtr space) {
  hid_t id;
  H5_CONSTRUCT(id, H5Acreate,
               (parent, name.c_str(), type->GetId(), space->GetId(),
                H5P_DEFAULT, H5P_DEFAULT));
  return std::make_shared<Attribute>(id);
}

AttributeSharedPtr Attribute::Open(hid_t parent, const std::string& name) {
  hid_t id;
  H5_CONSTRUCT(id, H5Aopen, (parent, name.c_str(), H5P_DEFAULT));
  return std::make_shared<Attribute>(id);
}

Attribute::~Attribute() {
  H5_CLOSETHIS(H5Aclose, true);
}
void Attribute::Close() {
  H5_CLOSETHIS(H5Aclose, false);
}

DataSpaceSharedPtr Attribute::GetSpace() const {
  return std::make_shared<DataSpace>(H5Aget_space(m_Id));
}

File::File(hid_t id) : Object(id) {}
FileSharedPtr File::Create(const std::string& filename,
                           unsigned mode,
                           PListSharedPtr createPL,
                           PListSharedPtr accessPL) {
  hid_t id;
  H5_CONSTRUCT(id, H5Fcreate,
               (filename.c_str(), mode, createPL->GetId(), accessPL->GetId()));
  return std::make_shared<File>(id);
}
FileSharedPtr File::Open(const std::string& filename,
                         unsigned mode,
                         PListSharedPtr accessPL) {
  hid_t id;
  H5_CONSTRUCT(id, H5Fopen, (filename.c_str(), mode, accessPL->GetId()));
  return std::make_shared<File>(id);
}

File::~File() {
  H5_CLOSETHIS(H5Fclose, true);
}

void File::Close() {
  H5_CLOSETHIS(H5Fclose, false);
}

hsize_t File::GetNumElements() {
  GroupSharedPtr root = OpenGroup("/");
  return root->GetNumElements();
}
Group::Group(hid_t id) : Object(id) {}

Group::~Group() {
  H5_CLOSETHIS(H5Gclose, true);
}

void Group::Close() {
  H5_CLOSETHIS(H5Gclose, false);
}

hsize_t Group::GetNumElements() {
  H5G_info_t info;
  H5_CALL(H5Gget_info, (m_Id, &info));
  return info.nlinks;
}

DataSet::DataSet(hid_t id) : Object(id) {}

DataSet::~DataSet() {
  H5_CLOSETHIS(H5Dclose, true);
}

void DataSet::Close() {
  H5_CLOSETHIS(H5Dclose, false);
}

DataSpaceSharedPtr DataSet::GetSpace() const {
  return std::make_shared<DataSpace>(H5Dget_space(m_Id));
}

}  // namespace hemelb::H5
