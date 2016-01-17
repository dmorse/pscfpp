#ifndef PSCF_BLOCK_STUB_H
#define PSCF_BLOCK_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
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
   };

} 
#endif
