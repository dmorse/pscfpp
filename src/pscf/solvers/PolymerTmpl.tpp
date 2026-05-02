#ifndef PSCF_POLYMER_TMPL_TPP
#define PSCF_POLYMER_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerTmpl.h"                   
#include <pscf/math/arithmetic.h>
#include <util/containers/Pair.h>

namespace Pscf {

   class Edge;
   using namespace Util;

   /*
   * Constructor.
   */
   template <class BT, class PT, typename WT>
   PolymerTmpl<BT,PT,WT>::PolymerTmpl()
    : PolymerSpeciesT(),
      blocks_()
   {  ParamComposite::setClassName("PolymerTmpl"); }

   /*
   * Allocate blocks array.
   */
   template <class BT, class PT, typename WT>
   void PolymerTmpl<BT,PT,WT>::allocateBlocks()
   {  blocks_.allocate(PolymerSpeciesT::nBlock()); }

   /*
   * Read blocks array from parameter file.
   */
   template <class BT, class PT, typename WT>
   void PolymerTmpl<BT,PT,WT>::readBlocks(std::istream& in)
   {
      const int nBlock = PolymerSpeciesT::nBlock();  
      ParamComposite::readDArray<BT>(in, "blocks", blocks_, nBlock); 
   }

   /*
   * Read parameter file block.
   */
   template <class BT, class PT, typename WT>
   void PolymerTmpl<BT,PT,WT>::readParameters(std::istream& in)
   {

      // Call PoymerSpecies base class member function
      // Initializes PolymerSpecies, and Edge and Vertex components
      PolymerSpeciesT::readParameters(in);

      // The remainder of this function sets and validates immutable
      // information about graph topology that is stored by propagators.

      PT * propagatorPtr = nullptr;
      PT const * sourcePtr = nullptr;
      Vertex const * headPtr = nullptr;
      Vertex const * tailPtr = nullptr;
      Pair<int> propId;
      int blockId, forwardId, reverseId, headId, tailId, i;
      const int nBlock = PolymerSpeciesT::nBlock();
      bool isHeadEnd, isTailEnd;

      // Set sources and end flags for all propagators
      for (blockId = 0; blockId < nBlock; ++blockId) {
         for (forwardId = 0; forwardId < 2; ++forwardId) {

            propagatorPtr = &block(blockId).propagator(forwardId);

            // Identify head and tail vertices
            if (forwardId == 0) {
              reverseId = 1;
            } else {
              reverseId = 0;
            }
            headId = block(blockId).vertexId(forwardId);
            tailId = block(blockId).vertexId(reverseId);
            headPtr = &PolymerSpeciesT::vertex(headId);  // head vertex
            tailPtr = &PolymerSpeciesT::vertex(tailId);  // tail vertex

            // Add pointers to source propagators
            for (i = 0; i < headPtr->size(); ++i) {
               propId = headPtr->inPropagatorId(i);
               if (propId[0] == blockId) {
                  UTIL_CHECK(propId[1] != forwardId);
               } else {
                  sourcePtr =
                     &block(propId[0]).propagator(propId[1]);
                  propagatorPtr->addSource(*sourcePtr);
               }
            }

            // Set vertex end flags
            isHeadEnd = (headPtr->size() == 1) ? true : false;
            isTailEnd = (tailPtr->size() == 1) ? true : false;
            propagatorPtr->setEndFlags(isHeadEnd, isTailEnd);

         } // end loop over forwardId (propagator direction id)
      } // end loop over blockId 

      // Check validity - throw Exception if error detected
      isValid();
   }

   /*
   * Checks validity of propagator data set in readParameters.
   *
   * This function only checks validity of propagator source and end flag
   * member data that is set in PolymerTmpl<BT,PT,WT>::readParameters. It 
   * does not check validity of members of PolymerSpecies, Edge, and Vertex
   * that are set and validated within the PolymerSpecies::readParameters 
   * base class member function. 
   */
   template <class BT, class PT, typename WT>
   void PolymerTmpl<BT,PT,WT>::isValid()
   {
      Vertex const * v0Ptr = nullptr;
      Vertex const * v1Ptr = nullptr;
      PT const * p0Ptr = nullptr;
      PT const * p1Ptr = nullptr;
      int bId, v0Id, v1Id;
      const int nVertex = PolymerSpeciesT::nVertex();
      const int nBlock = PolymerSpeciesT::nBlock();

      // Loop over blocks
      for (bId = 0; bId < nBlock; ++bId) {
         v0Id = block(bId).vertexId(0);
         v1Id = block(bId).vertexId(1);
         UTIL_CHECK(v0Id >= 0 && v0Id < nVertex);
         UTIL_CHECK(v1Id >= 0 && v1Id < nVertex);
         UTIL_CHECK(v0Id != v1Id);
         v0Ptr = &PolymerSpeciesT::vertex(v0Id);
         v1Ptr = &PolymerSpeciesT::vertex(v1Id);
         p0Ptr = &(block(bId).propagator(0));
         p1Ptr = &(block(bId).propagator(1));
         UTIL_CHECK(v0Ptr->size() >= 1);
         UTIL_CHECK(v1Ptr->size() >= 1);
         UTIL_CHECK(p0Ptr->nSource() == (v0Ptr->size() - 1));
         UTIL_CHECK(p1Ptr->nSource() == (v1Ptr->size() - 1));
         UTIL_CHECK(p0Ptr->isHeadEnd() == p1Ptr->isTailEnd());
         UTIL_CHECK(p0Ptr->isTailEnd() == p1Ptr->isHeadEnd());
         UTIL_CHECK(p0Ptr->isHeadEnd() == (v0Ptr->size() == 1));
         UTIL_CHECK(p1Ptr->isHeadEnd() == (v1Ptr->size() == 1));
      }

   }

   /*
   * Solve the MDE for all blocks of this polymer.
   */
   template <class BT, class PT, typename WT>
   void PolymerTmpl<BT,PT,WT>::solve(double phiTot)
   {
      const int nPropagator = PolymerSpeciesT::nPropagator();

      // Clear all propagators
      for (int j = 0; j < nPropagator; ++j) {
         propagator(j).setIsSolved(false);
      }

      // Solve modified diffusion equation for all propagators in
      // the order specified by function makePlan.
      for (int j = 0; j < nPropagator; ++j) {
         UTIL_CHECK(propagator(j).isReady());
         propagator(j).solve();
      }

      // Compute molecular partition function Q
      WT Q;
      block(0).propagator(0).computeQ(Q);

      // The PT::computeQ function returns an overall spatial average.
      // Correct for possible partial occupation of the unit cell.
      divEq(Q, phiTot);

      // Set q and compute phi or mu, depending on the ensemble
      SpeciesT::setQ(Q);

   }

   /*
   * Get a propagator, indexed by block and direction ids (non-const).
   */
   template <class BT, class PT, typename WT>
   PT& PolymerTmpl<BT,PT,WT>::propagator(int blockId, int directionId)
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a propagator, indexed by block and direction ids (const).
   */
   template <class BT, class PT, typename WT>
   PT const & 
   PolymerTmpl<BT,PT,WT>::propagator(int blockId, int directionId) const
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a propagator, indexed in order of computation.
   */
   template <class BT, class PT, typename WT>
   PT& PolymerTmpl<BT,PT,WT>::propagator(int id)
   {
      Pair<int> propId = PolymerSpeciesT::propagatorId(id);
      return propagator(propId[0], propId[1]);
   }

}
#endif
