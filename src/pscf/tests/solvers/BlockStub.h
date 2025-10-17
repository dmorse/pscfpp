#ifndef PSCF_BLOCK_STUB_H
#define PSCF_BLOCK_STUB_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PropagatorStub.h"
#include <pscf/solvers/BlockTmpl.tpp>
#include <util/containers/DArray.h>

namespace Pscf
{ 

   class BlockStub : public BlockTmpl<PropagatorStub, DArray<double> >
   {

   public:

      typedef PropagatorStub PropagatorT;
      typedef PropagatorStub::FieldT FieldT; 

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
      void setupSolver(FieldT const & wField)
      {}
   
      /**
      * Compute monomer concentration for this block.
      */ 
      void computeConcentration(double prefactor)
      {}
   
   };

}
#endif
