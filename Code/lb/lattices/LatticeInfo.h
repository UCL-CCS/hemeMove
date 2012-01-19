#ifndef HEMELB_LB_LATTICES_LATTICEINFO_H
#define HEMELB_LB_LATTICES_LATTICEINFO_H

#include "util/Vector3D.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      class LatticeInfo
      {
        public:
          inline LatticeInfo(unsigned numberOfVectors,
                             const util::Vector3D<int>* vectors,
                             const Direction* inverseVectorIndicesIn) :
              numVectors(numberOfVectors), vectorSet(), inverseVectorIndices()
          {
            for (Direction direction = 0; direction < numberOfVectors; ++direction)
            {
              vectorSet.push_back(util::Vector3D<int>(vectors[direction]));
              inverseVectorIndices.push_back(inverseVectorIndicesIn[direction]);
            }
          }

          inline unsigned GetNumVectors() const
          {
            return numVectors;
          }

          inline const util::Vector3D<int>& GetVector(unsigned index) const
          {
            return vectorSet[index];
          }

          inline unsigned GetInverseIndex(unsigned index) const
          {
            return inverseVectorIndices[index];
          }

        private:
          const unsigned numVectors;
          std::vector<util::Vector3D<int> > vectorSet;
          std::vector<Direction> inverseVectorIndices;
      };
    }
  }
}
#endif
