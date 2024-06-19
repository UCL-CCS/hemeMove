// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_OCT_H5_HPP
#define HLBGMYTOOL_OCT_H5_HPP

namespace hemelb::H5 {

template <class T>
typename DataTypeConversionPolicy<T>::ConvertedType
DataTypeConversionPolicy<T>::Convert(const T& obj) {
  return obj;
}

template <class T>
T DataTypeConversionPolicy<T>::Deconvert(
    const typename DataTypeConversionPolicy<T>::ConvertedType& obj) {
  return obj;
}

template <class T>
DataTypeSharedPtr DataTypeTraits<T>::GetType() {
  return PredefinedDataType::Native<T>();
}

template <class T>
const void* DataTypeTraits<T>::GetAddress(
    const DataTypeTraits<T>::ConvertedType& obj) {
  return &obj;
}

template <class T>
void* DataTypeTraits<T>::GetAddress(DataTypeTraits<T>::ConvertedType& obj) {
  return &obj;
}

template <class T>
typename DataTypeTraits<T>::ConvertedType DataTypeTraits<T>::Convert(
    const T& obj) {
  return Converter::Convert(obj);
}
template <class T>
T DataTypeTraits<T>::Deconvert(
    const typename DataTypeTraits<T>::ConvertedType& obj) {
  return Converter::Deconvert(obj);
}
template <>
struct DataTypeConversionPolicy<std::string> {
  static const bool MustConvert = true;
  typedef const char* ConvertedType;
  typedef const char* ConvertedVectorElemType;
  inline static ConvertedType Convert(const std::string& obj) {
    return obj.c_str();
  }
  inline static std::string Deconvert(const ConvertedType& obj) {
    return std::string(obj);
  }
};

template <>
inline DataTypeSharedPtr DataTypeTraits<std::string>::GetType() {
  return DataType::String();
}

template <class T>
DataTypeSharedPtr PredefinedDataType::Native() {
  return std::make_shared<PredefinedDataType>(DataTypeTraits<T>::NativeType);
}

template <typename T>
void CanHaveAttributes::SetAttribute(const std::string& name, const T& value) {
  DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
  DataSpaceSharedPtr space = DataSpace::Scalar();
  AttributeSharedPtr attr = CreateAttribute(name, type, space);

  typename DataTypeTraits<T>::ConvertedType conv =
      DataTypeTraits<T>::Convert(value);
  H5_CALL(H5Awrite,
          (attr->GetId(), type->GetId(), DataTypeTraits<T>::GetAddress(conv)));
}

template <typename T>
void CanHaveAttributes::SetAttribute(const std::string& name,
                                     const std::vector<T>& value) {
  typedef std::vector<
      typename DataTypeConversionPolicy<T>::ConvertedVectorElemType>
      Vec;
  Vec converted_vals;
  DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
  DataSpaceSharedPtr space = DataSpace::OneD(value.size());
  AttributeSharedPtr attr = CreateAttribute(name, type, space);

  const void* converted_buf = NULL;
  if (DataTypeConversionPolicy<T>::MustConvert) {
    converted_vals.resize(value.size());
    for (size_t i = 0; i < value.size(); ++i)
      converted_vals[i] = DataTypeConversionPolicy<T>::Convert(value[i]);
    converted_buf = &converted_vals[0];
  } else {
    converted_buf = &value[0];
  }

  H5_CALL(H5Awrite, (attr->GetId(), type->GetId(), converted_buf));
}

template <typename T>
void CanHaveAttributes::GetAttribute(const std::string& name, T& value) {
  DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
  DataSpaceSharedPtr space = DataSpace::Scalar();
  AttributeSharedPtr attr = OpenAttribute(name);

  typename DataTypeTraits<T>::ConvertedType conv;

  H5_CALL(H5Aread,
          (attr->GetId(), type->GetId(), DataTypeTraits<T>::GetAddress(conv)));
  value = DataTypeTraits<T>::Deconvert(conv);
}

template <typename T>
void CanHaveAttributes::GetAttribute(const std::string& name,
                                     std::vector<T>& value) {
  typedef std::vector<
      typename DataTypeConversionPolicy<T>::ConvertedVectorElemType>
      Vec;
  Vec converted_vals;
  DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
  AttributeSharedPtr attr = OpenAttribute(name);
  DataSpaceSharedPtr space = attr->GetSpace();
  if (H5Sget_simple_extent_ndims(space->GetId()) != 1)
    throw Error("vector data not 1D");
  hsize_t len, maxdim;
  H5Sget_simple_extent_dims(space->GetId(), &len, &maxdim);

  value.resize(len);
  void* converted_buf = NULL;
  if (DataTypeConversionPolicy<T>::MustConvert) {
    converted_vals.resize(len);
    converted_buf = &converted_vals[0];
  } else {
    converted_buf = &value[0];
  }

  H5_CALL(H5Aread, (attr->GetId(), type->GetId(), converted_buf));

  if (DataTypeConversionPolicy<T>::MustConvert) {
    typename Vec::iterator src = converted_vals.begin(),
                           end = converted_vals.end();
    typename std::vector<T>::iterator dest = value.begin();
    for (; src != end; ++src, ++dest) {
      *dest = DataTypeTraits<T>::Deconvert(*src);
    }
  }
}

template <class T>
DataSetSharedPtr CanHaveGroupsDataSets::CreateWriteDataSet(
    const std::string& name,
    const std::vector<T>& data,
    PListSharedPtr createPL,
    PListSharedPtr accessPL) {
  DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
  DataSpaceSharedPtr space = DataSpace::OneD(data.size());
  DataSetSharedPtr dataset =
      CreateDataSet(name, type, space, createPL, accessPL);
  dataset->Write(data);
  return dataset;
}

}  // namespace hemelb::H5

#endif
