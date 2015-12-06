#ifndef PFTS_SYSTEM_STUB_H
#define PFTS_SYSTEM_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/solvers/SystemTmpl.h>
#include <pfts/solvers/PolymerStub.h>
#include <pfts/solvers/SolventStub.h>

namespace Pfts{ 

   class SystemStub;

   class SystemStub 
    : public SystemTmpl<PolymerStub, SolventStub, DArray<double>, DArray<double> >
   {

   public:

      SystemStub();
   
   };

} 
#endif
