#ifndef PSCF_POLYMER_STUB_H
#define PSCF_POLYMER_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockStub.h"
#include <pscf/PolymerTmpl.h>

namespace Pscf
{ 

   class PolymerStub : public PolymerTmpl<BlockStub>
   {

   public:

      PolymerStub()
      {  setClassName("Polymer"); }

   };

} 
#endif
