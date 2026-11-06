#ifndef PSCF_POLYMER_STUB_H
#define PSCF_POLYMER_STUB_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockStub.h"
#include "PropagatorStub.h"
#include <pscf/solvers/PolymerTmpl.tpp>

namespace Pscf
{ 

   class PolymerStub : public PolymerTmpl<BlockStub, PropagatorStub, double>
   {

   public:

      PolymerStub()
      {  setClassName("Polymer"); }

   };

} 
#endif
