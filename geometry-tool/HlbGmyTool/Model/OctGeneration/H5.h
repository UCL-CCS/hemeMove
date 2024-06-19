// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_OCT_H5_H
#define HLBGMYTOOL_OCT_H5_H

#include <exception>
#include <memory>
#include <string>
#include <vector>

#include <hdf5.h>

namespace hemelb::H5 {
class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define H5_ERRORMSG(func) \
  "HDF5 error in API function " #func " " __FILE__ ":" TOSTRING(__LINE__)

#define H5_CONSTRUCT(ans, func, args) \
  {                                   \
    hid_t ret = func args;            \
    if (ret < 0)                      \
      throw Error(H5_ERRORMSG(func)); \
    ans = ret;                        \
  }

#define H5_CALL(func, args)           \
  {                                   \
    herr_t ret = func args;           \
    if (ret < 0)                      \
      throw Error(H5_ERRORMSG(func)); \
  }

#define H5_CALLNOTHROW(func, args)         \
  {                                        \
    herr_t ret = func args;                \
    if (ret < 0) {                         \
      std::cerr << H5_ERRORMSG(func) "\n"; \
      std::terminate();                    \
    }                                      \
  }

// Forward declare
class Object;
typedef std::shared_ptr<Object> ObjectSharedPtr;
class DataType;
typedef std::shared_ptr<DataType> DataTypeSharedPtr;
class DataSpace;
typedef std::shared_ptr<DataSpace> DataSpaceSharedPtr;
class CanHaveAttributes;
typedef std::shared_ptr<CanHaveAttributes> CanHaveAttributesSharedPtr;
class Attribute;
typedef std::shared_ptr<Attribute> AttributeSharedPtr;
class CanHaveGroupsDataSets;
typedef std::shared_ptr<CanHaveGroupsDataSets> CanHaveGroupsDataSetsSharedPtr;
class Group;
typedef std::shared_ptr<Group> GroupSharedPtr;
class File;
typedef std::shared_ptr<File> FileSharedPtr;
class DataSet;
typedef std::shared_ptr<DataSet> DataSetSharedPtr;
class PList;
typedef std::shared_ptr<PList> PListSharedPtr;

/// HDF5 base class
class Object {
 public:
  // Objects aren't copyable.
  Object(Object const&) = delete;
  Object& operator=(Object const&) = delete;
  virtual ~Object() noexcept;

  virtual void Close() = 0;
  inline hid_t GetId() const { return m_Id; }
  // Overload cast to the HDF5 ID type to make objects
  // transparently usable in the HDF5 API.
  inline operator hid_t() const { return GetId(); }

 protected:
  Object();
  explicit Object(hid_t id);
  hid_t m_Id;
};

// PropertyList objects
class PList : public Object {
 public:
  /// Default options
  static PListSharedPtr Default();
  /// Properties for object creation
  static PListSharedPtr ObjectCreate();
  /// Properties for file creation
  static PListSharedPtr FileCreate();
  /// Properties for file access
  static PListSharedPtr FileAccess();
  /// Properties for dataset creation
  static PListSharedPtr DatasetCreate();
  /// Properties for dataset access
  static PListSharedPtr DatasetAccess();
  /// Properties for raw data transfer
  static PListSharedPtr DatasetXfer();
  /// Properties for file mounting
  static PListSharedPtr FileMount();
  /// Properties for group creation
  static PListSharedPtr GroupCreate();
  /// Properties for group access
  static PListSharedPtr GroupAccess();
  /// Properties for datatype creation
  static PListSharedPtr DatatypeCreate();
  /// Properties for datatype access
  static PListSharedPtr DatatypeAccess();
  /// Properties for character encoding when encoding strings or object names
  static PListSharedPtr StringCreate();
  /// Properties for attribute creation
  static PListSharedPtr AttributeCreate();
  /// Properties governing the object copying process
  static PListSharedPtr ObjectCopy();
  /// Properties governing link creation
  static PListSharedPtr LinkCreate();
  /// Properties governing link traversal when accessing objects
  static PListSharedPtr LinkAccess();

  PList();
  ~PList() noexcept override;
  void Close() override;
  void SetChunk(const std::vector<hsize_t>& dims);
  void SetDeflate(const unsigned level = 1);
  // void SetMpio(Comm comm);
  // void SetDxMpioCollective();
  // void SetDxMpioIndependent();

  PList(hid_t cls);
};

// Input iterator of names of things (Links and Attributes).
//
// Need to supply a Policy class that provides types, behaviour and
// constants.
template <typename Policy>
class NameIterator {
 private:
  using ObjT = Policy::ObjectType;

  ObjT* m_obj;

  hsize_t m_idx;
  hsize_t m_next;
  hsize_t m_size;
  std::string m_currentName;

  static herr_t helper(hid_t id,
                       char const* name,
                       Policy::InfoT const* info,
                       void* op_data) {
    auto iter = static_cast<NameIterator*>(op_data);
    iter->m_currentName = name;
    return 1;
  }

  void next() {
    m_idx = m_next;
    if (m_idx < m_size) {
      H5_CALL(Policy::ITER_FUNC,
              (m_obj->GetId(), Policy::ITER_INDEX, Policy::ITER_ORDER, &m_next,
               NameIterator::helper, this));
    }
  }

 public:
  // Fulfill std::forward_iterator concept.
  using difference_type = hssize_t;
  using value_type = std::string const;
  using reference = std::string const&;
  using pointer = std::string const*;
  using iterator_category = std::input_iterator_tag;

  NameIterator(ObjT* obj, hsize_t idx)
      : m_obj{obj},
        /*m_idx{idx}, - initialised by next() below */
        m_next{idx},
        m_size{Policy::GetObjectSize(obj)} {
    next();
  }

  // Iterator gives the name of the current thing.
  const std::string& operator*() const { return m_currentName; }

  friend inline bool operator==(const NameIterator& lhs,
                                const NameIterator& rhs) {
    if (lhs.m_obj == rhs.m_obj) {
      if (lhs.m_idx == rhs.m_idx) {
        return true;
      }
    }
    return false;
  }

  friend inline bool operator!=(const NameIterator& lhs,
                                const NameIterator& rhs) {
    return !(lhs == rhs);
  }

  hsize_t GetPos() const { return m_idx; }

  NameIterator& operator++() {
    next();
    return *this;
  }
  NameIterator operator++(int) {
    NameIterator ans = *this;
    next();
    return ans;
  }
};

// Forward iterator over the names of links of a Group/File. These
// can be other Groups or DataSets.
//
// Iteration order is fastest available rather than sorted in any
// particular way.
struct LinkIterPolicy {
  using ObjectType = CanHaveGroupsDataSets;
  using InfoT = H5L_info_t;
  typedef herr_t (*IterFuncT)(hid_t,
                              H5_index_t,
                              H5_iter_order_t,
                              hsize_t*,
                              H5L_iterate_t,
                              void*);
  static constexpr IterFuncT ITER_FUNC = H5Literate;
  static constexpr auto ITER_INDEX = H5_INDEX_NAME;
  static constexpr auto ITER_ORDER = H5_ITER_NATIVE;
  // Definition below, once ObjectType is defined.
  static inline hsize_t GetObjectSize(ObjectType* o);
};
using LinkIterator = NameIterator<LinkIterPolicy>;

/// Mixin for objects that contain groups and datasets (Group and File)
class CanHaveGroupsDataSets : public virtual Object {
 public:
  // Create a group with the given name. The createPL can be
  // omitted to use the default properties.
  GroupSharedPtr CreateGroup(const std::string& name,
                             PListSharedPtr createPL = PList::Default(),
                             PListSharedPtr accessPL = PList::Default());

  // Create a dataset with the name, type and space.
  // The createPL can be omitted to use the defaults.
  DataSetSharedPtr CreateDataSet(const std::string& name,
                                 DataTypeSharedPtr type,
                                 DataSpaceSharedPtr space,
                                 PListSharedPtr createPL = PList::Default(),
                                 PListSharedPtr accessPL = PList::Default());

  // Create a dataset containing the data supplied
  // The createPL can be omitted to use the defaults
  template <class T>
  DataSetSharedPtr CreateWriteDataSet(
      const std::string& name,
      const std::vector<T>& data,
      PListSharedPtr createPL = PList::Default(),
      PListSharedPtr accessPL = PList::Default());

  // Open an existing group.
  // The accessPL can be omitted to use the defaults
  GroupSharedPtr OpenGroup(const std::string& name,
                           PListSharedPtr accessPL = PList::Default()) const;

  // Open an existing dataset
  // The accessPL can be omitted to use the defaults
  DataSetSharedPtr OpenDataSet(
      const std::string& name,
      PListSharedPtr accessPL = PList::Default()) const;

  // Get the number of links within this Group or File.
  virtual hsize_t GetNumElements() = 0;
  // Get iterator over the contained elements
  LinkIterator begin();
  LinkIterator end();
};

inline hsize_t LinkIterPolicy::GetObjectSize(ObjectType* o) {
  return o->GetNumElements();
}

// Forward iterator over the attribute names.
// Iteration order is whatever's fastest.
struct AttrIterPolicy {
  using ObjectType = CanHaveAttributes const;
  using InfoT = H5A_info_t;
  typedef herr_t (*IterFuncT)(hid_t,
                              H5_index_t,
                              H5_iter_order_t,
                              hsize_t*,
                              H5A_operator2_t,
                              void*);
  static inline IterFuncT ITER_FUNC = H5Aiterate2;
  static constexpr auto ITER_INDEX = H5_INDEX_CRT_ORDER;
  static constexpr auto ITER_ORDER = H5_ITER_NATIVE;
  // Definition below, once ObjectType is defined.
  static inline hsize_t GetObjectSize(ObjectType* o);
};
using AttrIterator = NameIterator<AttrIterPolicy>;

/// Mixin for objects that can have attributes (Group, DataSet, DataType)
class CanHaveAttributes : public virtual Object {
 public:
  AttributeSharedPtr CreateAttribute(const std::string& name,
                                     DataTypeSharedPtr type,
                                     DataSpaceSharedPtr space);
  AttributeSharedPtr OpenAttribute(const std::string& name);

  template <class T>
  void SetAttribute(const std::string& name, const T& value);
  template <class T>
  void SetAttribute(const std::string& name, const std::vector<T>& value);

  template <class T>
  void GetAttribute(const std::string& name, T& value);
  template <class T>
  void GetAttribute(const std::string& name, std::vector<T>& value);

  hsize_t GetNumAttr() const;
  AttrIterator attr_begin() const;
  AttrIterator attr_end() const;
};

inline hsize_t AttrIterPolicy::GetObjectSize(ObjectType* o) {
  return o->GetNumAttr();
}

/// HDF5 DataSpace wrapper
class DataSpace : public Object {
 public:
  static DataSpaceSharedPtr Null();
  static DataSpaceSharedPtr Scalar();
  static DataSpaceSharedPtr OneD(hsize_t size);

  DataSpace();
  DataSpace(hid_t id);
  DataSpace(hsize_t size, hsize_t max = H5S_UNLIMITED - 1);
  DataSpace(const std::vector<hsize_t>& dims);
  DataSpace(const std::vector<hsize_t>& dims,
            const std::vector<hsize_t>& max_dims);
  ~DataSpace() noexcept override;

  void Close() override;
  void SelectRange(const hsize_t start, const hsize_t count);
};

// Policy class for the DataTypesTraits controlling whether data
// has to be converted in anyway before writing. (It could perhaps
// be DataTypeTraitsTraits but that's too horrible a name.)
//
// Default policy is to not convert at all.
template <class T>
struct DataTypeConversionPolicy {
  static const bool MustConvert = false;
  typedef T ConvertedType;
  typedef T ConvertedVectorElemType;
  static ConvertedType Convert(const T& obj);
  static T Deconvert(const ConvertedType& obj);
};

/// Traits class for HDF5 data types.
template <class T>
struct DataTypeTraits {
  typedef DataTypeConversionPolicy<T> Converter;

  /***
   * Define this for a specialision for any HDF5 NATIVE type you want to use.
   * See http://hdfgroup.org/HDF5/doc/UG/UG_frame11Datatypes.html
   */
  static const hid_t NativeType;

  typedef typename Converter::ConvertedType ConvertedType;

  static ConvertedType Convert(const T& obj);
  static T Deconvert(const ConvertedType& obj);
  /**
   * Get the address of the start of the data.
   * Default implementation just uses "&"
   */
  static const void* GetAddress(const ConvertedType& obj);
  static void* GetAddress(ConvertedType& obj);
  /**
   * Return a DataType object representing T.
   * Default implementation just calls PredefinedDataType::Native<T>()
   */
  static DataTypeSharedPtr GetType();
  static DataTypeSharedPtr GetType(const T& obj);
};

/// Wrap and HDF5 data type object. Technically this can have attributes, but
/// not really bothered.
class DataType : public Object {
 public:
  static DataTypeSharedPtr String(size_t len = 0);
  template <class T>
  static DataTypeSharedPtr OfObject(const T& obj) {
    return DataTypeTraits<T>::GetType();
  }

  template <class T>
  static DataTypeSharedPtr Array(std::initializer_list<hsize_t> dims) {
    hsize_t nd = dims.size();
    hid_t ans;
    H5_CONSTRUCT(ans, H5Tarray_create,
                 (DataTypeTraits<T>::GetType()->GetId(), nd, dims.begin()));

    return std::make_shared<DataType>(ans);
  }

  DataType(hid_t id);
  ~DataType() noexcept override;
  void Close() override;
  DataTypeSharedPtr Copy() const;
};

/// Predefined HDF data types that must not be closed when done with.
class PredefinedDataType : public DataType {
 public:
  template <class T>
  static DataTypeSharedPtr Native();
  static DataTypeSharedPtr CS1();
  PredefinedDataType(hid_t);
  ~PredefinedDataType() noexcept override;
  void Close() override;
};

/// HDF5 Attribute Wrapper
class Attribute : public Object {
  // This type can be used in a range-for expression to iterate
  // attribute names - get one using the `Iterate` static member
  // function below.
  struct AttrIterationHelper {
    CanHaveAttributes const* thing;
    inline AttrIterator begin() const { return thing->attr_begin(); }
    inline AttrIterator end() const { return thing->attr_end(); }
  };

 public:
  Attribute(hid_t id);
  ~Attribute() noexcept override;
  void Close() override;
  DataSpaceSharedPtr GetSpace() const;
  inline static AttrIterationHelper Iterate(CanHaveAttributes const& thing) {
    return AttrIterationHelper{&thing};
  }

 private:
  static AttributeSharedPtr Create(hid_t parent,
                                   const std::string& name,
                                   DataTypeSharedPtr type,
                                   DataSpaceSharedPtr space);
  static AttributeSharedPtr Open(hid_t parent, const std::string& name);
  friend class CanHaveAttributes;
};

/// HDF5 file wrapper
class File : public CanHaveGroupsDataSets {
 public:
  static FileSharedPtr Create(const std::string& filename,
                              unsigned mode,
                              PListSharedPtr createPL = PList::Default(),
                              PListSharedPtr accessPL = PList::Default());
  static FileSharedPtr Open(const std::string& filename,
                            unsigned mode,
                            PListSharedPtr accessPL = PList::Default());
  ~File() noexcept override;
  void Close() override;
  hsize_t GetNumElements() override;

  File(hid_t id);
};

/// HDF5 Group wrapper
class Group : public CanHaveAttributes, public CanHaveGroupsDataSets {
 public:
  Group(hid_t id);
  ~Group() noexcept override;
  void Close() override;
  hsize_t GetNumElements() override;
  CanHaveAttributesSharedPtr operator[](hsize_t idx);
  CanHaveAttributesSharedPtr operator[](const std::string& key);
};

class DataSet : public CanHaveAttributes {
 public:
  DataSet(hid_t id);
  ~DataSet() noexcept override;
  void Close() override;
  DataSpaceSharedPtr GetSpace() const;

  template <class T>
  void Write(const std::vector<T>& data) {
    DataTypeSharedPtr mem_t = DataTypeTraits<T>::GetType();
    H5_CALL(H5Dwrite,
            (m_Id, mem_t->GetId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]));
  }
  template <class T>
  void Write(const std::vector<T>& data,
             DataSpaceSharedPtr filespace,
             PListSharedPtr dxpl = PList::Default()) {
    DataTypeSharedPtr mem_t = DataTypeTraits<T>::GetType();
    DataSpaceSharedPtr memspace = DataSpace::OneD(data.size());

    H5Dwrite(m_Id, mem_t->GetId(), memspace->GetId(), filespace->GetId(),
             dxpl->GetId(), &data[0]);
  }
  template <class T>
  void Read(std::vector<T>& data) {
    DataTypeSharedPtr mem_t = DataTypeTraits<T>::GetType();
    DataSpaceSharedPtr space = GetSpace();
    if (H5Sget_simple_extent_ndims(space->GetId()) != 1)
      throw Error("vector data not 1D");
    hsize_t len, maxdim;
    H5Sget_simple_extent_dims(space->GetId(), &len, &maxdim);

    data.resize(len);

    H5_CALL(H5Dread,
            (m_Id, mem_t->GetId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]));
  }
};

}  // namespace hemelb::H5

#include "H5.hpp"
#endif
