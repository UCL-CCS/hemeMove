// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_JUNKYANG_H
#define HEMELB_LB_STREAMERS_JUNKYANG_H

#include <set>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "hassert.h"
#include "units.h"
#include "lb/streamers/Common.h"

namespace hemelb::lb
{


    /**
       * Template to produce Streamers that can cope with fluid-fluid, fluid-solid
       * (using the Junk&Yang method, see below) and fluid-iolet links. Requires a Link Streamer
       * class that will handle the iolet links (or NullLink if wall-only).
       *
       * This class implements the Junk & Yang no-slip boundary condition as described in
       * M. Junk and Z. Yang "One-point boundary condition for the lattice Boltzmann method", Phys Rev E 72 (2005)
       */
    template<link_streamer IoletLinkImpl>
    class JunkYangFactory
    {
    public:
        using CollisionType = typename IoletLinkImpl::CollisionType;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename CollisionType::LatticeType;
        static constexpr bool has_iolet = !std::same_as<IoletLinkImpl, NullLink<CollisionType>>;
    private:
        using matrix = boost::numeric::ublas::matrix<distribn_t>;
        using identity_matrix = boost::numeric::ublas::identity_matrix<distribn_t>;
        using permutation_matrix = boost::numeric::ublas::permutation_matrix<std::size_t>;
        using c_vector = boost::numeric::ublas::c_vector<distribn_t, LatticeType::NUMVECTORS>;
        using vector = boost::numeric::ublas::vector<distribn_t>;
    public:
        JunkYangFactory(InitParams& initParams) :
              collider(initParams), bulkLinkDelegate(collider, initParams),
                  ioletLinkDelegate(collider, initParams), THETA(0.7)/*,
                  domainData(*initParams.latDat)*/
        {
            auto& dom = *initParams.latDat;
            for (auto&& [site_begin, site_end]: initParams.siteRanges)
            {
              for (site_t siteIdx = site_begin; siteIdx < site_end; ++siteIdx)
              {
                geometry::Site<const geometry::Domain> localSite = dom.GetSite(siteIdx);
//                    domainData.GetSite(siteIdx);
                // Only consider walls - the initParams .siteRanges should take care of that for us, but check anyway
                if (localSite.IsWall())
                {
                  ConstructVelocitySets(siteIdx, dom);
                  AssembleKMatrix(siteIdx, dom);
                  AssembleLMatrix(siteIdx);
                }
              }
            }
            FactoriseLMatrices();
          }

        ~JunkYangFactory()
        {
            // Free dynamically allocated permutation matrices.
            for (auto& luPermutationMatrix : luPermutationMatrices)
            {
                delete luPermutationMatrix.second;
            }
        }

        void StreamAndCollide(const site_t firstIndex, const site_t lastIndex,
                              const LbmParameters* lbmParams,
                              geometry::FieldData& latticeData,
                              lb::MacroscopicPropertyCache& propertyCache)
        {
            for (site_t siteIndex = firstIndex; siteIndex < lastIndex; siteIndex++)
            {
              HASSERT(latticeData.GetSite(siteIndex).IsWall());
              HASSERT(lMatrices.find(siteIndex) != lMatrices.end());

              auto&& site = latticeData.GetSite(siteIndex);

              VarsType hydroVars(site);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              // In the first step, we stream and collide as we would for the BulkStreamer
              // streamer.
              collider.CalculatePreCollision(hydroVars, site);
              collider.Collide(lbmParams, hydroVars);

              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                if (site.HasWall(direction))
                {
                  // Do nothing for now.
                }
                else if (site.HasIolet(direction))
                {
                    if constexpr (has_iolet) {
                        ioletLinkDelegate.StreamLink(lbmParams, latticeData, site, hydroVars, direction);
                    } else {
                        throw (Exception() << "No iolet type configured but iolet link encountered");
                    }
                }
                else
                {
                  bulkLinkDelegate.StreamLink(lbmParams, latticeData, site, hydroVars, direction);
                }
              }

              // The following for loops prepare the data required by DoPostStep

              // TODO: Ideally, this should be done in the loop over directions above
              // but the f's are permuted by the below making it tricky.
              unsigned index = 0;
              for (auto incomingVelocityIter = incomingVelocities[siteIndex].begin();
                  incomingVelocityIter != incomingVelocities[siteIndex].end();
                  ++incomingVelocityIter, ++index)
              {
                int inverseDirection = LatticeType::INVERSEDIRECTIONS[*incomingVelocityIter];

                fPostCollision[siteIndex](index) =
                    hydroVars.GetFPostCollision()[*incomingVelocityIter];
                fPostCollisionInverseDir[siteIndex](index) =
                    hydroVars.GetFPostCollision()[inverseDirection];
                fOld[siteIndex](index) = site.GetFOld<LatticeType>()[*incomingVelocityIter];
              }

              for (auto outgoingVelocityIter = outgoingVelocities[siteIndex].begin();
                  outgoingVelocityIter != outgoingVelocities[siteIndex].end();
                  ++outgoingVelocityIter, ++index)
              {
                fPostCollision[siteIndex](index) =
                    hydroVars.GetFPostCollision()[*outgoingVelocityIter];
                fOld[siteIndex](index) = site.GetFOld<LatticeType>()[*outgoingVelocityIter];
              }

                UpdateCachePostCollision(site,
                                         hydroVars,
                                         lbmParams,
                                         propertyCache);
            }

        }

        void PostStep(const site_t firstIndex, const site_t lastIndex,
                      const LbmParameters* lbmParams, geometry::FieldData& latticeData,
                      lb::MacroscopicPropertyCache& propertyCache)
        {
            for (site_t siteIndex = firstIndex; siteIndex < lastIndex; siteIndex++)
            {
              HASSERT(latticeData.GetSite(siteIndex).IsWall());
              HASSERT(lMatrices.find(siteIndex) != lMatrices.end());
              HASSERT(luPermutationMatrices.find(siteIndex) != luPermutationMatrices.end());

              vector rVector;
              AssembleRVector(siteIndex, latticeData,rVector);

              // assemble RHS
              vector systemRHS = fPostCollisionInverseDir[siteIndex] - rVector;

              // lu_substitue will overwrite the RHS with the solution
              vector &systemSolution = systemRHS;
              lu_substitute(lMatrices[siteIndex],
                            *luPermutationMatrices[siteIndex],
                            systemSolution);

              // Update the distribution function for incoming velocities with the solution of the linear system
              unsigned index = 0;
              for (auto incomingVelocityIter = incomingVelocities[siteIndex].begin();
                  incomingVelocityIter != incomingVelocities[siteIndex].end();
                  ++incomingVelocityIter, ++index)
              {
                * (latticeData.GetFNew(siteIndex * LatticeType::NUMVECTORS + *incomingVelocityIter)) =
                    systemSolution[index];
              }

              auto&& site = latticeData.GetSite(siteIndex);
              for (unsigned int outgoingVelocityIter : outgoingVelocities[siteIndex])
              {
                  if constexpr (has_iolet) {
                      if (site.HasIolet(outgoingVelocityIter)) {
                          ioletLinkDelegate.PostStepLink(latticeData, site, outgoingVelocityIter);
                      }
                  }
              }
            }
        }

    private:
          CollisionType collider;
          BulkLink<CollisionType> bulkLinkDelegate;
          IoletLinkImpl ioletLinkDelegate;
          //! Problem dimension (2D, 3D)
          static const unsigned DIMENSION = 3U;
          //! Vector coordinate arbitrarily chosen in the paper
          static const unsigned ALPHA = DIMENSION - 1;
          //! theta constant in the theta-method used for interpolation (0 for fully explicit, 1 for fully implicit)
          const distribn_t THETA;

          //! Reference to the lattice object used for initialisation
          //const geometry::Domain& domainData;

          //! Sets of incoming velocities (those with an inverse direction crossing a wall boundary)
          std::map<site_t, std::set<Direction> > incomingVelocities;
          //! Sets of outgoing velocities (the complement of incomingVelocities)
          std::map<site_t, std::set<Direction> > outgoingVelocities;

          //! Map containing the K matrices used to assemble both the left-hand- and the right-hand-side of the linear systems for each site
          std::map<site_t, matrix> kMatrices;
          //! Map containing the linear system left-hand-sides for each site
          std::map<site_t, matrix > lMatrices;
          /**
           * Map containing the permutations generated by the LU factorisation of each linear system. The permutaation matrices are dinamically allocated
           * because there isn't a default constructor for ublas::permutation_matrix (unlike other ublas data structures) and the size is not known until
           * run time.
           */
          std::map<site_t, permutation_matrix*> luPermutationMatrices;

          //! Map of postcollision distributions with incoming/outgoing ordering.
          std::map<site_t, c_vector > fPostCollision;
          //! Map of postcollision distributions for the inverse directions of those in the incoming set of velocities.
          std::map<site_t, vector > fPostCollisionInverseDir;
          //! Map of distributions in the previous time step with incoming/outgoing ordering.
          std::map<site_t, c_vector > fOld;

          /**
           * Construct the incoming/outgoing velocity sets for site siteLocalIndex
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           */
          inline void ConstructVelocitySets(site_t contiguousSiteIndex, geometry::Domain const& dom)
          {
            auto site = dom.GetSite(contiguousSiteIndex);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              int inverseDirection = LatticeType::INVERSEDIRECTIONS[direction];
              if (site.HasWall(inverseDirection))
              {
                incomingVelocities[contiguousSiteIndex].insert(direction);
              }
              else
              {
                outgoingVelocities[contiguousSiteIndex].insert(direction);
              }
            }
            fPostCollisionInverseDir[contiguousSiteIndex].resize(incomingVelocities[contiguousSiteIndex].size());
          }

          /**
           * Assemble the K matrix for site contiguousSiteIndex.
           *
           * K is a rectangular matrix (num_incoming_vels x LatticeType::NUMVECTORS). Our
           * implementation places the columns corresponding to the set of incoming velocities
           * first followed by those corresponding to outgoing velocities. This makes easier
           * handling the subblocks with uBLAS.
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           */
          inline void AssembleKMatrix(site_t contiguousSiteIndex, geometry::Domain const& dom)
          {
            geometry::Site<const geometry::Domain> site =
                dom.GetSite(contiguousSiteIndex);

            unsigned incomingVelsSetSize = incomingVelocities[contiguousSiteIndex].size();

            kMatrices[contiguousSiteIndex].resize(incomingVelsSetSize, LatticeType::NUMVECTORS);

            /*
             *  Assemble K(:, 0:num_incoming_velocities-1)
             */
            unsigned rowIndex = 0;
            for (auto rowIndexIncomingVelocity = incomingVelocities[contiguousSiteIndex].begin();
                rowIndexIncomingVelocity != incomingVelocities[contiguousSiteIndex].end();
                ++rowIndexIncomingVelocity, ++rowIndex)
            {
              // |c_i|^2, where c_i is the i-th velocity vector
              const int rowRowdirectionsInnProd = LatticeType::CX[*rowIndexIncomingVelocity]
                  * LatticeType::CX[*rowIndexIncomingVelocity]
                  + LatticeType::CY[*rowIndexIncomingVelocity]
                      * LatticeType::CY[*rowIndexIncomingVelocity]
                  + LatticeType::CZ[*rowIndexIncomingVelocity]
                      * LatticeType::CZ[*rowIndexIncomingVelocity];

              unsigned columnIndex = 0;
              for (auto columnIndexIncomingVelocity = incomingVelocities[contiguousSiteIndex].begin();
                  columnIndexIncomingVelocity != incomingVelocities[contiguousSiteIndex].end();
                  ++columnIndexIncomingVelocity, ++columnIndex)
              {
                // |c_i|^2, where c_i is the i-th velocity vector
                const int colColdirectionsInnProd = LatticeType::CX[*columnIndexIncomingVelocity]
                    * LatticeType::CX[*columnIndexIncomingVelocity]
                    + LatticeType::CY[*columnIndexIncomingVelocity]
                        * LatticeType::CY[*columnIndexIncomingVelocity]
                    + LatticeType::CZ[*columnIndexIncomingVelocity]
                        * LatticeType::CZ[*columnIndexIncomingVelocity];

                // (c_i \dot c_j)^2, where c_{i,j} are the {i,j}-th velocity vectors
                const int rowColdirectionsInnProd = LatticeType::CX[*rowIndexIncomingVelocity]
                    * LatticeType::CX[*columnIndexIncomingVelocity]
                    + LatticeType::CY[*rowIndexIncomingVelocity]
                        * LatticeType::CY[*columnIndexIncomingVelocity]
                    + LatticeType::CZ[*rowIndexIncomingVelocity]
                        * LatticeType::CZ[*columnIndexIncomingVelocity];

                HASSERT(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction ) *rowIndexIncomingVelocity])
                    >= 0);
                HASSERT(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction ) *rowIndexIncomingVelocity])
                    < 1);

                kMatrices[contiguousSiteIndex](rowIndex, columnIndex) =
                    (-3.0 / 2.0)
                        * (3.0
                            - 6
                                * site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity]))
                        * LatticeType::EQMWEIGHTS[*rowIndexIncomingVelocity]
                        * ( (rowColdirectionsInnProd * rowColdirectionsInnProd)
                            - (rowRowdirectionsInnProd / 3.0)
                            - LatticeType::VECTORS[*rowIndexIncomingVelocity][ALPHA]
                                * LatticeType::VECTORS[*rowIndexIncomingVelocity][ALPHA]
                                * (colColdirectionsInnProd - (DIMENSION / 3.0)));

                HASSERT(fabs(kMatrices[contiguousSiteIndex](rowIndex, columnIndex)) < 1e3);
              }

              /*
               * Assemble K(:, num_incoming_velocities:num_incoming_velocities+num_outgoing_velocities-1)
               */
              for (auto columnIndexOutgoingVelocity =
                  outgoingVelocities[contiguousSiteIndex].begin();
                  columnIndexOutgoingVelocity != outgoingVelocities[contiguousSiteIndex].end();
                  ++columnIndexOutgoingVelocity, ++columnIndex)
              {
                // |c_i|^2, where c_i is the i-th velocity vector
                const int colColdirectionsInnProd = LatticeType::CX[*columnIndexOutgoingVelocity]
                    * LatticeType::CX[*columnIndexOutgoingVelocity]
                    + LatticeType::CY[*columnIndexOutgoingVelocity]
                        * LatticeType::CY[*columnIndexOutgoingVelocity]
                    + LatticeType::CZ[*columnIndexOutgoingVelocity]
                        * LatticeType::CZ[*columnIndexOutgoingVelocity];

                // (c_i \dot c_j)^2, where c_{i,j} are the {i,j}-th velocity vectors
                const int rowColdirectionsInnProd = LatticeType::CX[*rowIndexIncomingVelocity]
                    * LatticeType::CX[*columnIndexOutgoingVelocity]
                    + LatticeType::CY[*rowIndexIncomingVelocity]
                        * LatticeType::CY[*columnIndexOutgoingVelocity]
                    + LatticeType::CZ[*rowIndexIncomingVelocity]
                        * LatticeType::CZ[*columnIndexOutgoingVelocity];

                HASSERT(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction ) *rowIndexIncomingVelocity])
                    >= 0);
                HASSERT(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction ) *rowIndexIncomingVelocity])
                    < 1);

                kMatrices[contiguousSiteIndex](rowIndex, columnIndex) =
                    (-3.0 / 2.0)
                        * (3.0
                            - 6
                                * site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity]))
                        * LatticeType::EQMWEIGHTS[*rowIndexIncomingVelocity]
                        * ( (rowColdirectionsInnProd * rowColdirectionsInnProd)
                            - (rowRowdirectionsInnProd / 3.0)
                            - LatticeType::VECTORS[*rowIndexIncomingVelocity][ALPHA]
                                * LatticeType::VECTORS[*rowIndexIncomingVelocity][ALPHA]
                                * (colColdirectionsInnProd - (DIMENSION / 3.0)));

                HASSERT(fabs(kMatrices[contiguousSiteIndex](rowIndex, columnIndex)) < 1e3);
              }

            }
          }

          /**
           * Assemble the L matrix for site contiguousSiteIndex. L is a square matrix (num_incoming_vels x num_incoming_vels)
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           */
          inline void AssembleLMatrix(site_t siteLocalIndex)
          {
            unsigned incomingVelsSetSize = incomingVelocities[siteLocalIndex].size();

            lMatrices[siteLocalIndex] = identity_matrix(incomingVelsSetSize)
                + THETA
                    * subrange(kMatrices[siteLocalIndex],
                                      0,
                                      incomingVelsSetSize,
                                      0,
                                      incomingVelsSetSize);
          }

          /**
           * Compute the LU factorisation of all the L matrices
           */
          inline void FactoriseLMatrices()
          {
            for (auto iter = lMatrices.begin(); iter != lMatrices.end(); ++iter)
            {
              luPermutationMatrices[iter->first] =
                  new permutation_matrix(incomingVelocities[iter->first].size());
              // If this assertion trips, lMatrices[siteLocalIndex] is singular.
              HASSERT(lu_factorize(iter->second, *luPermutationMatrices[iter->first]) == 0);
            }
          }

          /**
           * Assemble the sigma vector required to assemble the system RHS. We are not including
           * the forcing term used in the paper to drive the flow. This might become necessary
           * for biocolloids.
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           * @param sigmaVector sigma vector
           */
          inline void AssembleSigmaVector(
              const site_t contiguousSiteIndex,
              c_vector& sigmaVector) // const
          {
            HASSERT(fPostCollision.find(contiguousSiteIndex) != fPostCollision.end());
            HASSERT(fOld.find(contiguousSiteIndex) != fOld.end());
            sigmaVector = fPostCollision[contiguousSiteIndex]
                - (1 - THETA) * fOld[contiguousSiteIndex];
          }

          /**
           * Assemble the r vector required to assemble the system RHS
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           * @param rVector r vector
           */
          inline void AssembleRVector(const site_t contiguousSiteIndex,
                                      geometry::FieldData const& fieldData,
                                      vector& rVector) //const
          {
            c_vector sigmaVector;
            AssembleSigmaVector(contiguousSiteIndex, sigmaVector);

            unsigned incomingVelsSetSize = incomingVelocities[contiguousSiteIndex].size();
            unsigned outgoingVelsSetSize = outgoingVelocities[contiguousSiteIndex].size();

            vector fNew(outgoingVelsSetSize);

            /*
             *  Assemble a vector with the updated values of the distribution function for the outgoing velocities,
             *  which have already been streamed
             */
            HASSERT(fNew.size() == outgoingVelocities[contiguousSiteIndex].size());
            unsigned index = 0;
            for (auto outgoingDirIter =
                outgoingVelocities[contiguousSiteIndex].begin();
                outgoingDirIter != outgoingVelocities[contiguousSiteIndex].end();
                ++outgoingDirIter, ++index)
            {
              fNew[index] = *fieldData.GetFNew(contiguousSiteIndex * LatticeType::NUMVECTORS
                  + *outgoingDirIter);
            }

            rVector = THETA
                * prod(subrange(kMatrices[contiguousSiteIndex],
                                              0,
                                              incomingVelsSetSize,
                                              incomingVelsSetSize,
                                              LatticeType::NUMVECTORS),
                              fNew) + prod(kMatrices[contiguousSiteIndex], sigmaVector);
          }
      };
}
#endif
