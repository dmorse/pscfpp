#ifndef PSCF_BLOCK_STUB_H
#define PSCF_BLOCK_STUB_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PropagatorStub.h"
#include <pscf/BlockTmpl.h>

namespace Pscf
{ 

   class BlockStub : public BlockTmpl<PropagatorStub>
   {

   public:

      typedef PropagatorStub Propagator;
      typedef PropagatorStub::WField WField; 

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
      void setupSolver(WField const & wField)
      {}
   
      /**
      * Compute monomer concentration for this block.
      */ 
      void computeConcentration(double prefactor)
      {}
   
   };

}
#endif
