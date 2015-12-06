#ifndef PFTS_PROPAGATOR_STUB_H
#define PFTS_PROPAGATOR_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/solvers/PropagatorTmpl.h>

namespace Pfts{ 

   class PropagatorStub;

   class PropagatorStub : public PropagatorTmpl<PropagatorStub>
   {

   public:

      PropagatorStub(){}
   
   };

} 
#endif
