#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerSpecies.h>    // base class
#include <util/containers/DArray.h>      // member template

// Classes used in implementation
#include <pscf/chem/Vertex.h>           
#include <pscf/chem/PolymerType.h>       
#include <util/containers/Pair.h>       

#include <cmath>

namespace Pscf
{

   class Edge;
   using namespace Util;

   /**
   * Descriptor and MDE solver for an acyclic branched block polymer.
   *
   * A PolymerTmpl<Block> object has an array of Block objects and an
   * array of Vertex objects that it inherits from the PolymerSpecies
   * base class.  Each Block has two propagator MDE solver objects 
   * associated with the two directions along each block. The solve() 
   * member function solves the modified diffusion equation (MDE) for 
   * all propagators in the molecule (i.e., all blocks, in both 
   * directions) and computes monomer concentration fields for all 
   * blocks.
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
      * function q, phi and mu are all set, and concentration fields
      * associated with all blocks in the polymer have been computed.
      * The implementation of PolymerTmpl::solve() the calls the solve()
      * function for all propagators in the molecule in a predetermined
      * order.
      *
      * The concrete subclass of PolymerTmpl<Block> in each implemenation
      * level namespace, which is named Polymer by convention, defines a 
      * function named "compute" that calls PolymerTmpl<Block>::solve. 
      * This compute function takes an array of chemical potential 
      * fields (w-fields) as an argument. Before calling the solve() 
      * function declared here, the compute() function of each such Polymer
      * class must pass the w-fields and any other mutable required data 
      * to all Block objects in order to set up the solution of the MDE 
      * within each block. 
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
      * a Block as a reference to an Edge (a block descriptor).
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
      using PolymerSpecies::edge;
      using PolymerSpecies::vertex;
      using PolymerSpecies::propagatorId;
      using PolymerSpecies::path;
      using PolymerSpecies::nBlock;
      using PolymerSpecies::nVertex;
      using PolymerSpecies::nPropagator;
      using PolymerSpecies::length;
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

   /*
   * Get a specified Edge (block descriptor) by non-const reference.
   */
   template <class Block>
   inline Edge& PolymerTmpl<Block>::edge(int id)
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
   * Get a Propagator indexed by block and direction (non-const).
   */
   template <class Block>
   inline
   typename Block::Propagator&
   PolymerTmpl<Block>::propagator(int blockId, int directionId)
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a Propagator indexed by block and direction (const).
   */
   template <class Block>
   inline
   typename Block::Propagator const &
   PolymerTmpl<Block>::propagator(int blockId, int directionId) const
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a propagator indexed in order of computation.
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
      Pair<int> propagatorId;
      int blockId, directionId, vertexId, i;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         // Add sources
         for (directionId = 0; directionId < 2; ++directionId) {
            vertexId = block(blockId).vertexId(directionId);
            vertexPtr = &vertex(vertexId);
            propagatorPtr = &block(blockId).propagator(directionId);
            for (i = 0; i < vertexPtr->size(); ++i) {
               propagatorId = vertexPtr->inPropagatorId(i);
               if (propagatorId[0] == blockId) {
                  UTIL_CHECK(propagatorId[1] != directionId);
               } else {
                  sourcePtr =
                     &block(propagatorId[0]).propagator(propagatorId[1]);
                  propagatorPtr->addSource(*sourcePtr);
               }
            }
         }

      }

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

      // Compute block concentration fields
      double prefactor = phi_ / ( q_ * length() );
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }

   }

}
#endif
