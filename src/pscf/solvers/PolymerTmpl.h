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
   * A PolymerTmpl<Block> object has arrays of Block and Vertex objects.
   * Each Block has two propagator MDE solver objects.  The solve() 
   * member function solves the modified diffusion equation (MDE) for 
   * all propagators in molecule.
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
      * This function should be called within a member function named
      * "compute" of each concrete subclass.  This compute function must
      * take an array of chemical potential fields (w-fields) as an
      * argument. Before calling the solve() function declared here,
      * the compute() function of each such subclass must setup the 
      * solvers by storing information required to solve the MDE within 
      * each block within a set of implementation dependent private 
      * members of the class that represents a block of a block polymer. 
      * After calling this solve() function, the compute function must
      * loop over blocks to compute the monomer concentration fields 
      * for all blocks.
      *
      * The optional parameter phiTot is only relevant to problems
      * such as thin films in which the material is excluded from part
      * of the unit cell by an inhogeneous constraint on the sum of
      * monomer concentration (i.e., a "mask").
      *
      * \param phiTot  fraction of volume occupied by material
      */
      virtual void solve(double phiTot = 1.0);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      virtual Edge& edge(int id);

      /**
      * Get a specified Edge (block descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      virtual Edge const& edge(int id) const;

      /**
      * Get a specified Block.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Block& block(int id);

      /**
      * Get a specified Block by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Block const& block(int id) const ;

      /**
      * Get propagator for a specific block and direction (non-const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator& propagator(int blockId, int directionId);

      /**
      * Get a propagator for a specific block and direction (const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator const & propagator(int blockId, int directionId) const;

      /**
      * Get propagator indexed in order of computation (non-const).
      *
      * The propagator index must satisfy 0 <= id < 2*nBlock.
      *
      * \param id propagator index, in order of computation plan
      */
      Propagator& propagator(int id);

      ///@}

      using PolymerSpecies::edge;
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
      * Allocate array of blocks.
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
