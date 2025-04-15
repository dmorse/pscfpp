#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerSpecies.h>    // base class
#include <pscf/chem/Species.h>           // base class
#include <util/param/ParamComposite.h>   // base class

#include <util/containers/Pair.h>        // member template
#include <util/containers/DArray.h>      // member template

#include <pscf/chem/Edge.h>
#include <pscf/chem/Vertex.h>
#include <pscf/chem/PolymerType.h>
#include <pscf/chem/PolymerModel.h>

namespace Pscf
{

   using namespace Util;

   /**
   * Descriptor and MDE solver for a block polymer.
   *
   * A PolymerTmpl<Block> object has an array of Block objects and an
   * array of Vertex objects that it inherits from the PolymerSpecies base
   * class.  Each Block has two propagator MDE solver objects associated
   * with the two directions along each block of a block polymer.
   *
   * The solve() member function solves the modified diffusion equation
   * (MDE) for all propagators in the molecule (i.e., all blocks, in both
   * directions).
   *
   * Each implementation-level subclass of Pscf contains a class named
   * Block that derived from Pscf::Edge, and a class named Polymer that
   * is derived from PolymerTmpl<Block>.
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
      * Solve modified diffusion equation.
      *
      * Upon return, propagators for all blocks, molecular partition
      * function q, phi and mu are all set.  The implementation of
      * PolymerTmpl::solve() the calls the solve() function for all
      * propagators in the molecule in a predeternined order, then
      * computes q, and finally computes mu from phi or phi from mu,
      * depending on the ensemble.
      *
      * The concrete subclass of PolymerTmpl<Block> in each implemenation
      * level namespace, which is named Polymer by convention, defines 
      * a function named "compute" that calls PolymerTmpl<Block>::solve.
      * This compute function takes an array of chemical potential fields
      * (w-fields) as an argument. Before calling the solve() function
      * declared here, the compute() function of each such Polymer class
      * must pass the w-fields and any other required mutable data to 
      * all Block objects in order to set up the solution of the MDE 
      * within each block. After calling the solve() function, the 
      * Polymer::compute function must also compute concentrations for 
      * all blocks, each of which is owned by associated Block objects.
      *
      * The optional parameter phiTot is only relevant to problems
      * involving a mask that excludes material from part of the unit
      * cell, as for thin film problems.
      *
      * \param phiTot  fraction of unit cell volume occupied by material
      */
      virtual void solve(double phiTot = 1.0);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * The edge member functions implement pure virtual functions
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
      * Get a specified Block (block solver and descriptor)
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
      * Get propagator indexed in order of computation (non-const).
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
      virtual void allocateBlocks();

      /**
      * Read array of data for blocks from parameter file
      *
      * \param in parameter input stream
      */
      virtual void readBlocks(std::istream& in);

   private:

      /// Array of Block objects in this polymer.
      DArray<Block> blocks_;

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

   template <class Block>
   void PolymerTmpl<Block>::allocateBlocks()
   {  blocks_.allocate(nBlock()); }

   template <class Block>
   void PolymerTmpl<Block>::readBlocks(std::istream& in)
   {  readDArray<Block>(in, "blocks", blocks_, nBlock()); }

   template <class Block>
   void PolymerTmpl<Block>::readParameters(std::istream& in)
   {

      PolymerSpecies::readParameters(in);

      // Set sources for all propagators
      Vertex const * vertexPtr = 0;
      Propagator const * sourcePtr = 0;
      Propagator * propagatorPtr = 0;
      Pair<int> propId;
      int blockId, directionId, vertexId, i;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         // Add sources
         for (directionId = 0; directionId < 2; ++directionId) {
            vertexId = block(blockId).vertexId(directionId);
            vertexPtr = &vertex(vertexId);
            propagatorPtr = &block(blockId).propagator(directionId);
            for (i = 0; i < vertexPtr->size(); ++i) {
               propId = vertexPtr->inPropagatorId(i);
               if (propId[0] == blockId) {
                  UTIL_CHECK(propId[1] != directionId);
               } else {
                  sourcePtr =
                     &block(propId[0]).propagator(propId[1]);
                  propagatorPtr->addSource(*sourcePtr);
               }
            } // end loop over inputs to head vertex
         } // end loop over directionId
      } // end loop over blockId

      // If using bead model, set propagator vertex ownership
      Block* blockPtr;
      if (PolymerModel::isBead()) {
         bool own0, own1;
         for (int blockId = 0; blockId < nBlock(); ++blockId) {
            blockPtr = &(blocks_[blockId]);
            own0 = blockPtr->ownsVertex(0);
            own1 = blockPtr->ownsVertex(1);
            blockPtr->propagator(0).setVertexOwnership(own0, own1);
            blockPtr->propagator(1).setVertexOwnership(own1, own0);
         }
      }

      // If using bead model, check ownership of vertex beads
      if (PolymerModel::isBead()) {
         int vSize, ib, id, nOwner;
         Pair<int> propId;

         // Loop over vertices
         for (int vertexId = 0; vertexId < nVertex(); ++vertexId) {
            vSize = vertex(vertexId).size();

            // Incoming propagators
            nOwner = 0;
            for (int ip = 0; ip < vSize; ++ip) {
               propId = vertex(vertexId).inPropagatorId(ip);
               ib = propId[0];
               id = propId[1];
               if (block(ib).propagator(id).ownsTail()) {
                 ++nOwner;
               }
            }
            UTIL_CHECK(nOwner == 1);

            // Outgoing propagators
            nOwner = 0;
            for (int i = 0; i < vSize; ++i) {
               propId = vertex(vertexId).outPropagatorId(i);
               ib = propId[0];
               id = propId[1];
               if (block(ib).propagator(id).ownsHead()) {
                 ++nOwner;
               }
            }
            UTIL_CHECK(nOwner == 1);

         } // end loop over vertices

         // Check consistency of block and propagator ownership flags
         bool own0, own1;
         for (int blockId = 0; blockId < nBlock(); ++blockId) {
            blockPtr = &(block(blockId));
            own0 = blockPtr->ownsVertex(0);
            own1 = blockPtr->ownsVertex(1);
            UTIL_CHECK(blockPtr->propagator(0).ownsHead() == own0);
            UTIL_CHECK(blockPtr->propagator(1).ownsTail() == own0);
            UTIL_CHECK(blockPtr->propagator(0).ownsTail() == own1);
            UTIL_CHECK(blockPtr->propagator(1).ownsHead() == own1);
         }

      } // end if (PolymerModel::isBead())

   }

   /*
   * Compute a solution to the MDE and block concentrations.
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

      // Compute molecular partition function q_
      q_ = block(0).propagator(0).computeQ();

      // The Propagator::computeQ function returns a spatial average.
      // Correct for partial occupation of the unit cell.
      q_ = q_/phiTot;

      // Compute mu_ or phi_, depending on ensemble
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q_);
      } else
      if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q_;
      }

      #if 0
      // Compute block concentration fields
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi_ / ( q_ * length() );
      } else
      if (PolymerModel::isBead()) {
         prefactor = phi_ / ( q_ * (double)nBead() );
      }
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }
      #endif

   }

}
#endif
