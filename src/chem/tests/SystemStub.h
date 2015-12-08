#ifndef CHEM_SYSTEM_STUB_H
#define CHEM_SYSTEM_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerStub.h"
#include "SolventStub.h"
#include <chem/SystemTmpl.h>

namespace Chem{ 

   class SystemStub;

   class SystemStub 
    : public SystemTmpl<PolymerStub, SolventStub>
   {

   public:

      SystemStub()
      {  setClassName("System"); }
   
   };

} 
#endif
