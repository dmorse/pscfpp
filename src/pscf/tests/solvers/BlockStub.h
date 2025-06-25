#ifndef PSCF_BLOCK_STUB_H
#define PSCF_BLOCK_STUB_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PropagatorStub.h"
#include <pscf/solvers/BlockTmpl.h>

namespace Pscf
{ 

   class BlockStub : public BlockTmpl<PropagatorStub>
   {

   public:

      typedef PropagatorStub PropagatorT;
      typedef PropagatorStub::WFieldT WFieldT; 

      /**
      * Constructor.
      */
      BlockStub()
      {
         propagator(0).setBlock(*this);
         propagator(1).setBlock(*this);
      }

      /**
      * Compute monomer concentration for this block.
      */ 
      void setupSolver(WFieldT const & wField)
      {}
   
      /**
      * Compute monomer concentration for this block.
      */ 
      void computeConcentration(double prefactor)
      {}
   
   };

}
#endif
