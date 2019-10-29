#ifndef PSCF_SOLVENT_STUB_H
#define PSCF_SOLVENT_STUB_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/SolventTmpl.h>
#include "PropagatorStub.h"

namespace Pscf
{ 

   class SolventStub : public SolventTmpl<PropagatorStub>
   {

   public:

      SolventStub(){}

      void compute(PropagatorStub::WField const &)
      {}

   };

} 
#endif
