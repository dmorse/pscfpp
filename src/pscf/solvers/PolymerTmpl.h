#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerSpecies.h>    // base class
#include <util/param/ParamComposite.h>   // base class

#include <util/containers/Pair.h>        // member template
#include <util/containers/DArray.h>      // member template

#include <pscf/chem/Vertex.h>
#include <pscf/chem/PolymerType.h>
#include <pscf/chem/PolymerModel.h>

namespace Pscf
{

   class Edge;
   using namespace Util;

   /**
   * Template for an MDE solver and descriptor for a block polymer.
   *
   * A PolymerTmpl<Block> object has an array of Block objects, as well
   * as an array of Vertex objects inherited from the PolymerSpecies base
   * class.  Each Block has two Propagator MDE solver objects associated
   * with the two directions along each block.
   *
   * The solve() member function solves the modified diffusion equation
   * (MDE) for all propagators in the molecule (i.e., all blocks, in both
   * directions) in a pre-defined order.
   *
   * Each implementation-level sub-namespace of Pscf contains a class
   * named Block that is derived from Pscf::Edge, and a class named
   * Polymer that is derived from PolymerTmpl<Block>.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class Block>
   class PolymerTmpl : public PolymerSpecies
   {

   public:

      // Modified diffusion equation solver for one block.
      typedef typename Block::Propagator Propagator;

      // Monomer concentration field.
      typedef typename Propagator::CField CField;

      // Chemical potential field.
      typedef typename Propagator::WField WField;

      /**
      * Constructor.
      */
      PolymerTmpl();

      /**
      * Destructor.
      */
      ~PolymerTmpl();

      /**
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Solve modified diffusion equation for all propagators.
      *
      * Upon return, all propagators for all blocks, molecular partition
      * function q, phi and mu are all set.  The implementation of
      * PolymerTmpl::solve() the calls the solve() function for all
      * propagators in the molecule in a predeternined order, then
      * computes q, and finally computes mu from phi or phi from mu,
      * depending on the ensemble (open or closed) for this species.
      * This solve() function does not compute monomer concentration
      * field.
      *
      * Each implementation-level namespace defines a concrete subclass
      * of PolymerTmpl<Block> that is named Polymer by convention. Each
      * such Polymer class defines a function named "compute" that takes
      * an array of chemical fields (w-fields) as an argument, and that
      * calls PolymerTmpl<Block>::solve internally.  Before calling the
      * solve() function declared here, the Polymer::compute() function
      * must pass the w-fields and any other required mutable data to all
      * Block objects in order to set up the MDE solver for each block.
      * After calling the solve() function, the Polymer::compute function
      * must also compute monomer concentrations for all blocks, each of
      * which is stored in a field container owned by the associated Block.
      *
      * The optional parameter phiTot is only relevant to problems
      * involving a Mask that excludes material from part of the unit
      * cell, as done to define thin film problems.
      *
      * \param phiTot  fraction of unit cell volume occupied by material
      */
      virtual void solve(double phiTot = 1.0);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * The edge member functions implements pure virtual functions
      * defined by the PolymerSpecies base class, and provide access to
      * each Block as a reference to an Edge (a block descriptor).
      *
      * \param id block index, 0 <= id < nBlock
      */
      Edge& edge(int id) final;

      /**
      * Get a specified Edge (block descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Edge const& edge(int id) const final;

      /**
      * Get a specified Block (block solver and descriptor).
      *
      * \param id block index, 0 <= id < nBlock
      */
      Block& block(int id);

      /**
      * Get a specified Block (solver and descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Block const& block(int id) const ;

      /**
      * Get the Propagator for a specific block and direction (non-const).
      *
      * For an edge that terminates at vertices with vertex indices given
      * by the return values of Edge::vertexId(0) and Edge::vertexId(1):
      *
      *    - direction 0 propagates from vertexId(0) to vertexId(1)
      *    - direction 1 propagates from vertexId(1) to vertexId(0)
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator& propagator(int blockId, int directionId);

      /**
      * Get the Propagator for a specific block and direction (const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator const & propagator(int blockId, int directionId) const;

      /**
      * Get a propagator indexed in order of computation (non-const).
      *
      * The propagator index argument must satisfy 0 <= id < 2*nBlock.
      *
      * \param id  propagator index, in order of computation plan
      */
      Propagator& propagator(int id);

      ///@}

      // Inherited public members
      using PolymerSpecies::vertex;
      using PolymerSpecies::propagatorId;
      using PolymerSpecies::path;
      using PolymerSpecies::nBlock;
      using PolymerSpecies::nVertex;
      using PolymerSpecies::nPropagator;
      using PolymerSpecies::length;
      using PolymerSpecies::nBead;
      using PolymerSpecies::type;

      using Species::phi;
      using Species::mu;
      using Species::q;
      using Species::ensemble;

   protected:

      /**
      * Allocate array of Block objects.
      */
      void allocateBlocks() final;

      /**
      * Read array of data for blocks from parameter file
      *
      * \param in parameter input stream
      */
      void readBlocks(std::istream& in) final;

   private:

      /// Array of Block objects in this polymer.
      DArray<Block> blocks_;

      /**
      * Check validity of internal data structures.
      *
      * An Exception is thrown if any error is detected.
      */
      void isValid();

   };

   // Inline functions

   /*
   * Get a specified Edge (block descriptor) by non-const reference.
   */
   template <class Block>
   inline
   Edge& PolymerTmpl<Block>::edge(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Edge (block descriptor) by const reference.
   */
   template <class Block>
   inline
   Edge const & PolymerTmpl<Block>::edge(int id) const
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by non-const reference.
   */
   template <class Block>
   inline
   Block& PolymerTmpl<Block>::block(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by const reference.
   */
   template <class Block>
   inline
   Block const & PolymerTmpl<Block>::block(int id) const
   {  return blocks_[id]; }

   /*
   * Get a Propagator, indexed by block and direction ids (non-const).
   */
   template <class Block>
   inline
   typename Block::Propagator&
   PolymerTmpl<Block>::propagator(int blockId, int directionId)
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a Propagator, indexed by block and direction ids (const).
   */
   template <class Block>
   inline
   typename Block::Propagator const &
   PolymerTmpl<Block>::propagator(int blockId, int directionId) const
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a Propagator, indexed in order of computation.
   */
   template <class Block>
   inline
   typename Block::Propagator& PolymerTmpl<Block>::propagator(int id)
   {
      Pair<int> propId = propagatorId(id);
      return propagator(propId[0], propId[1]);
   }

   // Non-inline functions

   /*
   * Constructor.
   */
   template <class Block>
   PolymerTmpl<Block>::PolymerTmpl()
    : PolymerSpecies(),
      blocks_()
   {  setClassName("PolymerTmpl"); }

   /*
   * Destructor.
   */
   template <class Block>
   PolymerTmpl<Block>::~PolymerTmpl()
   {}

   /*
   * Allocate blocks array.
   */
   template <class Block>
   void PolymerTmpl<Block>::allocateBlocks()
   {  blocks_.allocate(nBlock()); }

   /*
   * Read blocks array from parameter file.
   */
   template <class Block>
   void PolymerTmpl<Block>::readBlocks(std::istream& in)
   {  readDArray<Block>(in, "blocks", blocks_, nBlock()); }

   /*
   * Read parameter file block.
   */
   template <class Block>
   void PolymerTmpl<Block>::readParameters(std::istream& in)
   {

      // Call base class method
      PolymerSpecies::readParameters(in);

      Vertex const * headPtr = nullptr;
      Vertex const * tailPtr = nullptr;
      Propagator const * sourcePtr = nullptr;
      Propagator * propagatorPtr = nullptr;
      Pair<int> propId;
      int blockId, forwardId, reverseId, headId, tailId, i;
      bool isHeadEnd, isTailEnd;

      // Set sources and end flags for all propagators
      for (blockId = 0; blockId < nBlock(); ++blockId) {
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
            headPtr = &vertex(headId);
            tailPtr = &vertex(tailId);

            // Add sources
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
            isHeadEnd = false;
            isTailEnd = false;
            if (headPtr->size() == 1) isHeadEnd = true;
            if (tailPtr->size() == 1) isTailEnd = true;
            propagatorPtr->setEndFlags(isHeadEnd, isTailEnd);

         } // end loop over forwardId (propagator direction id)
      } // end loop over blockId 

      // Check validity - throw Exception if error detected
      isValid();
   }

   template <class Block>
   void PolymerTmpl<Block>::isValid()
   {
      Vertex const * v0Ptr = nullptr;
      Vertex const * v1Ptr = nullptr;
      Propagator const * p0Ptr = nullptr;
      Propagator const * p1Ptr = nullptr;
      int bId, v0Id, v1Id;

      // Check propagator context properties
      for (bId = 0; bId < nBlock(); ++bId) {
         v0Id = block(bId).vertexId(0);
         v1Id = block(bId).vertexId(1);
         UTIL_CHECK(v0Id >= 0 && v0Id < nVertex());
         UTIL_CHECK(v1Id >= 0 && v1Id < nVertex());
         UTIL_CHECK(v0Id != v1Id);
         v0Ptr = &vertex(v0Id);
         v1Ptr = &vertex(v1Id);
         p0Ptr = &(block(bId).propagator(0));
         p1Ptr = &(block(bId).propagator(1));
         UTIL_CHECK(p0Ptr->nSource() == (v0Ptr->size() - 1));
         UTIL_CHECK(p1Ptr->nSource() == (v1Ptr->size() - 1));
         UTIL_CHECK(p0Ptr->isHeadEnd() == p1Ptr->isTailEnd());
         UTIL_CHECK(p0Ptr->isHeadEnd() == (v0Ptr->size() == 1));
         UTIL_CHECK(p1Ptr->isHeadEnd() == (v1Ptr->size() == 1));
      }

   }

   /*
   * Solve the MDE for all blocks.
   */
   template <class Block>
   void PolymerTmpl<Block>::solve(double phiTot)
   {

      // Clear all propagators
      for (int j = 0; j < nPropagator(); ++j) {
         propagator(j).setIsSolved(false);
      }

      // Solve modified diffusion equation for all propagators in
      // the order specified by function makePlan.
      for (int j = 0; j < nPropagator(); ++j) {
         UTIL_CHECK(propagator(j).isReady());
         propagator(j).solve();
      }

      // Compute molecular partition function Q
      double Q = block(0).propagator(0).computeQ();

      // The Propagator::computeQ function returns a spatial average.
      // Correct for partial occupation of the unit cell.
      Q = Q/phiTot;

      // Set q and compute phi or mu, depending on the ensemble
      Species::setQ(Q);

   }

}
#endif
