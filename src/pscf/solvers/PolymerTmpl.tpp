#ifndef PSCF_POLYMER_TMPL_TPP
#define PSCF_POLYMER_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerTmpl.h"                   
#include <pscf/chem/PolymerModel.h>
#include <pscf/math/arithmetic.h>
#include <util/containers/Pair.h>

#include <cmath>
#include <limits>

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
   {
      assign(qShifted_, 0.0);
      ParamComposite::setClassName("PolymerTmpl");
   }

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

      // Compute partition function Q of the propagators, which are
      // solutions in the shifted fields w - wShift of all blocks
      WT Q;
      block(0).propagator(0).computeQ(Q);

      // The PT::computeQ function returns an overall spatial average.
      // Correct for possible partial occupation of the unit cell.
      divEq(Q, phiTot);
      assign(qShifted_, Q);

      // Remove the field shifts: the shift of a block multiplies Q by
      // exp(wShift*length) (thread) or exp(wShift*nBead) (bead), so
      // ln q = ln Q - sum over blocks of these exponents
      double lnShift = 0.0;
      for (int j = 0; j < PolymerSpeciesT::nBlock(); ++j) {
         if (PolymerModel::isThread()) {
            lnShift += block(j).wShift() * block(j).length();
         } else {
            lnShift += block(j).wShift() * (double) block(j).nBead();
         }
      }
      WT lnQ;
      assignLog(lnQ, Q);
      subEq(lnQ, lnShift);

      // Set ln q and compute phi or mu, depending on the ensemble
      SpeciesT::setLnQ(lnQ);

   }

   /*
   * Raise the field shifts of all blocks if qShifted() is too small.
   */
   template <class BT, class PT, typename WT>
   void PolymerTmpl<BT,PT,WT>::raiseShifts(
                   std::vector<double> & shifts,
                   std::function<void(int, double)> const & setupBlock,
                   std::function<double(int)> const & solveBlock,
                   double phiTot)
   {
      // A raise by logLimit in ln q cannot make a q below qMin overflow
      const double logLimit
                  = std::log(std::numeric_limits<double>::max()) - 32.0;
      const double qMin = std::exp(-logLimit);
      if (!(qShifted_ < qMin)) return;

      const int nb = PolymerSpeciesT::nBlock();
      const bool isThread = PolymerModel::isThread();

      // Raise shifts until an average q of shifted propagators is one,
      // given the initial q, the length l over which the shifts apply,
      // and a function raise(c) that adds c to the shifts concerned,
      // solves, and returns the new q
      auto raiseToOne = [&](double q, double l,
                            std::function<double(double)> const & raise) {
         for (int k = 0; k < 1000 && q < qMin; ++k) {
            q = raise(logLimit/l);
         }
         if (q >= qMin && std::isfinite(q)) {
            raise(-std::log(q)/l);
         }
      };

      // Calibrate each block alone
      const std::vector<double> initial(shifts);
      for (int j = 0; j < nb; ++j) {
         auto raiseBlock = [&](double c) -> double {
            shifts[j] += c;
            setupBlock(j, shifts[j]);
            return solveBlock(j);
         };
         const double l = isThread ? block(j).length()
                                   : (double) block(j).nBead();
         raiseToOne(raiseBlock(0.0), l, raiseBlock);
      }

      // Solve the polymer, and raise all shifts by a common constant
      const double n = isThread ? PolymerSpeciesT::length()
                                : (double) PolymerSpeciesT::nBead();
      auto raisePolymer = [&](double c) -> double {
         for (int j = 0; j < nb; ++j) {
            shifts[j] += c;
            setupBlock(j, shifts[j]);
         }
         solve(phiTot);
         return qShifted_;
      };
      raiseToOne(raisePolymer(0.0), n, raisePolymer);

      // If the result is not representable, raise the initial shifts
      // by a common constant only
      if (!(qShifted_ >= qMin && std::isfinite(qShifted_))) {
         shifts = initial;
         raiseToOne(raisePolymer(0.0), n, raisePolymer);
      }
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
