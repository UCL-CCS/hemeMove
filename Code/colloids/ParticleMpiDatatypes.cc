// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/Particle.h"
//#include "units.h"

namespace hemelb
{
  namespace net
  {
    // boiler-plate template specialisation for colloids Particle object
    template<>
    MPI_Datatype MpiDataTypeTraits<colloids::Particle>::RegisterMpiDataType()
    {
      return colloids::Particle().CreateMpiDatatypeWithVelocity();
    }

    // boiler-plate template specialisation for colloids PersistedParticle object
    template<>
    MPI_Datatype MpiDataTypeTraits<colloids::PersistedParticle>::RegisterMpiDataType()
    {
      return colloids::Particle().CreateMpiDatatypeWithPosition();
    }

    // boiler-plate template specialisation for incoming particle velocity
    template<>
    MPI_Datatype MpiDataTypeTraits<std::pair<unsigned long, util::Vector3D<double> > >::RegisterMpiDataType()
    {
      // MPI_Get_address specifies non-const pointers
      // so need to create non-const std::pair object
      std::pair<unsigned long, util::Vector3D<double> > temp;

      // we want a re-usable MPI data type
      // so we need relative displacements
      MPI_Aint baseAddress;
      MPI_Get_address(&temp, &baseAddress);

      // we have chosen to make each block of fields contain a single field
      // so, the number of field blocks is the same as the number of fields
      int numberOfFieldBlocks = 2;

      // and the length of every field block is one
      int lengthOfEachFieldBlock[] = { 1, 1 };

      // there is no guarantee that the fields will be contiguous, so
      // the displacement of each field must be determined separately
      MPI_Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
      MPI_Get_address(& (temp.first), & (displacementOfEachFieldBlock[0]));
      MPI_Get_address(& (temp.second), & (displacementOfEachFieldBlock[1]));
      displacementOfEachFieldBlock[0] -= baseAddress;
      displacementOfEachFieldBlock[1] -= baseAddress;

      // the built-in MPI datatype of each field must match the C++ type
      MPI_Datatype datatypeOfEachFieldBlock[] = { MPI_UNSIGNED_LONG, MpiDataType<
                                                      util::Vector3D<double> >() };

      // create a first draft of the MPI datatype for a Particle
      // the lower bound and displacements of fields are correct
      // but the extent may not include the whole derived object
      // specifically, we aren't sending velocity and bodyForces
      MPI_Datatype pairType;
      MPI_Type_create_struct(numberOfFieldBlocks,
                             lengthOfEachFieldBlock,
                             displacementOfEachFieldBlock,
                             datatypeOfEachFieldBlock,
                             &pairType);

      // obtain the current lower bound for the MPI datatype
      MPI_Aint lowerBound, extent;
      MPI_Type_get_extent(pairType, &lowerBound, &extent);

      // we can determine the actual extent of a Particle object
      // by concatenating two of them, using a contiguous vector
      // and finding the difference between their base addresses
      std::vector<std::pair<unsigned long, util::Vector3D<double> > > tempVectorOfPair(2, temp);
      MPI_Get_address(& (tempVectorOfPair[0]), &baseAddress);
      MPI_Get_address(& (tempVectorOfPair[1]), &extent);
      extent -= baseAddress;

      // resize the uncommitted first draft MPI datatype
      // with the current lower bound and the new extent
      MPI_Datatype pairType2;
      MPI_Type_create_resized(pairType, lowerBound, extent, &pairType2);

      // commit the MPI datatype and return it
      MPI_Type_commit(&pairType2);
      return pairType2;
    }
  }
  namespace colloids
  {

    /** creates a derived MPI datatype that represents a single particle object
     *  note - this data type uses displacements rather than absolute addresses
     *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
     *  when you no longer need this type, remember to call MPI_Type_free
     */
    MPI_Datatype Particle::CreateMpiDatatypeWithPosition() const
    {
      // MPI_Get_address specifies non-const pointers
      // so need a non-const copy of a particle object
      Particle temp(*this);

      // we want a re-usable MPI data type
      // so we need relative displacements
      MPI_Aint baseAddress;
      MPI_Get_address(&temp, &baseAddress);

      // we have chosen to make each block of fields contain a single field
      // so, the number of field blocks is the same as the number of fields
      int numberOfFieldBlocks = 8;

      // and the length of every field block is one
      int lengthOfEachFieldBlock[] = { 1, 1, 1, 1, 1, 1, 1, 1 };

      // there is no guarantee that the fields will be contiguous, so
      // the displacement of each field must be determined separately
      MPI_Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
      MPI_Get_address(& (temp.particleId), & (displacementOfEachFieldBlock[0]));
      MPI_Get_address(& (temp.ownerRank), & (displacementOfEachFieldBlock[1]));
      MPI_Get_address(& (temp.smallRadius_a0), & (displacementOfEachFieldBlock[2]));
      MPI_Get_address(& (temp.largeRadius_ah), & (displacementOfEachFieldBlock[3]));
      MPI_Get_address(& (temp.mass), & (displacementOfEachFieldBlock[4]));
      MPI_Get_address(& (temp.lastCheckpointTimestep), & (displacementOfEachFieldBlock[5]));
      MPI_Get_address(& (temp.markedForDeletionTimestep), & (displacementOfEachFieldBlock[6]));
      MPI_Get_address(& (temp.globalPosition), & (displacementOfEachFieldBlock[7]));
      for (int ndx = 0; ndx < numberOfFieldBlocks; ndx++)
        displacementOfEachFieldBlock[ndx] -= baseAddress;

      // the built-in MPI datatype of each field must match the C++ type
      MPI_Datatype datatypeOfEachFieldBlock[] = { MPI_UNSIGNED_LONG,
      MPI_INT,
                                                  MPI_DOUBLE,
                                                  MPI_DOUBLE,
                                                  MPI_DOUBLE,
                                                  MPI_UNSIGNED_LONG,
                                                  MPI_UNSIGNED_LONG, net::MpiDataType<
                                                      util::Vector3D<double> >() };

      // create a first draft of the MPI datatype for a Particle
      // the lower bound and displacements of fields are correct
      // but the extent may not include the whole derived object
      // specifically, we aren't sending velocity and bodyForces
      MPI_Datatype particleType;
      MPI_Type_create_struct(numberOfFieldBlocks,
                             lengthOfEachFieldBlock,
                             displacementOfEachFieldBlock,
                             datatypeOfEachFieldBlock,
                             &particleType);

      // obtain the current lower bound for the MPI datatype
      MPI_Aint lowerBound, extent;
      MPI_Type_get_extent(particleType, &lowerBound, &extent);

      // we can determine the actual extent of a Particle object
      // by concatenating two of them, using a contiguous vector
      // and finding the difference between their base addresses
      std::vector<Particle> tempVectorOfParticle(2, temp);
      MPI_Get_address(& (tempVectorOfParticle[0]), &baseAddress);
      MPI_Get_address(& (tempVectorOfParticle[1]), &extent);
      extent -= baseAddress;

      // resize the uncommitted first draft MPI datatype
      // with the current lower bound and the new extent
      MPI_Datatype particleType2;
      MPI_Type_create_resized(particleType, lowerBound, extent, &particleType2);

      // commit the MPI datatype and return it
      MPI_Type_commit(&particleType2);
      return particleType2;
    }

    /** creates a derived MPI datatype that represents a single particle object
     *  note - this data type uses displacements rather than absolute addresses
     *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
     *  when you no longer need this type, remember to call MPI_Type_free
     */
    MPI_Datatype Particle::CreateMpiDatatypeWithVelocity() const
    {
      // MPI_Get_address specifies non-const pointers
      // so need a non-const copy of a particle object
      Particle temp(*this);

      // we want a re-usable MPI data type
      // so we need relative displacements
      MPI_Aint baseAddress;
      MPI_Get_address(&temp, &baseAddress);

      // we have chosen to make each block of fields contain a single field
      // so, the number of field blocks is the same as the number of fields
      int numberOfFieldBlocks = 2;

      // and the length of every field block is one
      int lengthOfEachFieldBlock[] = { 1, 1 };

      // there is no guarantee that the fields will be contiguous, so
      // the displacement of each field must be determined separately
      MPI_Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
      MPI_Get_address(& (temp.particleId), & (displacementOfEachFieldBlock[0]));
      MPI_Get_address(& (temp.velocity), & (displacementOfEachFieldBlock[1]));
      displacementOfEachFieldBlock[0] -= baseAddress;
      displacementOfEachFieldBlock[1] -= baseAddress;

      // the built-in MPI datatype of each field must match the C++ type
      MPI_Datatype datatypeOfEachFieldBlock[] = { MPI_UNSIGNED_LONG, net::MpiDataType<
                                                      util::Vector3D<double> >() };

      // create a first draft of the MPI datatype for a Particle
      // the lower bound and displacements of fields are correct
      // but the extent may not include the whole derived object
      // specifically, we aren't sending velocity and bodyForces
      MPI_Datatype particleType;
      MPI_Type_create_struct(numberOfFieldBlocks,
                             lengthOfEachFieldBlock,
                             displacementOfEachFieldBlock,
                             datatypeOfEachFieldBlock,
                             &particleType);

      // obtain the current lower bound for the MPI datatype
      MPI_Aint lowerBound, extent;
      MPI_Type_get_extent(particleType, &lowerBound, &extent);

      // we can determine the actual extent of a Particle object
      // by concatenating two of them, using a contiguous vector
      // and finding the difference between their base addresses
      std::vector<Particle> tempVectorOfParticle(2, temp);
      MPI_Get_address(& (tempVectorOfParticle[0]), &baseAddress);
      MPI_Get_address(& (tempVectorOfParticle[1]), &extent);
      extent -= baseAddress;

      // resize the uncommitted first draft MPI datatype
      // with the current lower bound and the new extent
      MPI_Datatype particleType2;
      MPI_Type_create_resized(particleType, lowerBound, extent, &particleType2);

      // commit the MPI datatype and return it
      MPI_Type_commit(&particleType2);
      return particleType2;
    }

  }
}

