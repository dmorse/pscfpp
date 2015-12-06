#ifndef PFTS_POLYMER_STUB_H
#define PFTS_POLYMER_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/solvers/PolymerTmpl.h>
#include <pfts/solvers/PropagatorStub.h>

namespace Pfts{ 

   class PolymerStub;

   class PolymerStub : public PolymerTmpl<PropagatorStub, DArray<double> >
   {

   public:

      PolymerStub();
   
   };

} 
#endif
