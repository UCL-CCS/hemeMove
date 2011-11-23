#include "Vector3D.h"
#include "mpiInclude.h"

namespace hemelb
{
  template<typename vectorType>
  MPI_Datatype GenerateTypeForVector()
  {
    const int typeCount = 1;
    int blocklengths[typeCount] = { 3 };

    MPI_Datatype types[typeCount] = { MpiDataType<vectorType> () };

    MPI_Aint displacements[typeCount] = { 0 };

    MPI_Datatype ret;

    MPI_Type_struct(typeCount, blocklengths, displacements, types, &ret);

    MPI_Type_commit(&ret);
    return ret;
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<float> >::RegisterMpiDataType()
  {
    return GenerateTypeForVector<float> ();
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<site_t> >::RegisterMpiDataType()
  {
    return GenerateTypeForVector<site_t> ();
  }
}
