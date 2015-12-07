#ifndef PFTS_PROPAGATOR_STUB_H
#define PFTS_PROPAGATOR_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/PropagatorTmpl.h>
#include <util/containers/DArray.h>

namespace Pfts{ 

   class PropagatorStub;

   class PropagatorStub : public PropagatorTmpl<PropagatorStub>
   {

   public:

      PropagatorStub(){}

      void compute(const DArray<double>& wField)
      {
         setIsComplete(true);
      };
   
   };

} 
#endif
