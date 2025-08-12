#ifndef PSCF_PROPAGATOR_STUB_H
#define PSCF_PROPAGATOR_STUB_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h>
#include <util/containers/DArray.h>

namespace Pscf
{ 

   class PropagatorStub;
   class BlockStub;

   class PropagatorStub : public PropagatorTmpl<PropagatorStub>
   {

   public:

      // Public typedefs

      /**
      * Field type.
      */
      using FieldT = DArray<double>;

      // Member functions

      /**
      * Constructor.
      */
      PropagatorStub(){}

      /**
      * Associate this propagator with a block.
      *
      * \param block associated Block object.
      */ 
      void setBlock(BlockStub& block);

      /**
      * Get the associated Block object by reference.
      */
      BlockStub & block();

      /**
      * Solve the modified diffusion equation for this block.
      *
      * \param w chemical potential field for appropriate monomer type
      */
      void solve()
      {  setIsSolved(true); };
  
      /**
      * Compute monomer concentration for this block.
      */ 
      void computeConcentration(double prefactor)
      {};
   
      /**
      * Compute and return molecular partition function.
      */ 
      double computeQ()
      {  return 1.0; };

   private:

      /// Pointer to associated Block.
      BlockStub* blockPtr_;

   };

   /*
   * Associate this propagator with a block and direction
   */
   inline void PropagatorStub::setBlock(BlockStub& block)
   {  blockPtr_ = &block; }

   /*
   * Get the associated Block object.
   */
   inline BlockStub& PropagatorStub::block()
   {  return *blockPtr_; }

} 
#endif
