#ifndef PSCF_POLYMER_STUB_H
#define PSCF_POLYMER_STUB_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockStub.h"
#include "PropagatorStub.h"
#include <pscf/solvers/PolymerTmpl.tpp>

namespace Pscf
{ 

   class PolymerStub : public PolymerTmpl<BlockStub, PropagatorStub>
   {

   public:

      PolymerStub()
      {  setClassName("Polymer"); }

   };

} 
#endif
