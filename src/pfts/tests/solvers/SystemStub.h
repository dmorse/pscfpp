#ifndef PFTS_SYSTEM_STUB_H
#define PFTS_SYSTEM_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerStub.h"
#include "SolventStub.h"
#include <pfts/solvers/SystemTmpl.h>

namespace Pfts{ 

   class SystemStub;

   class SystemStub 
    : public SystemTmpl<PolymerStub, SolventStub, DArray<double>, DArray<double> >
   {

   public:

      SystemStub()
      {  setClassName("System"); }
   
   };

} 
#endif
