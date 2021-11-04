#ifndef PSCF_SOLVENT_STUB_H
#define PSCF_SOLVENT_STUB_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PropagatorStub.h"
#include <pscf/chem/Species.h>

namespace Pscf
{ 

   class SolventStub : public Species
   {

   public:

      SolventStub(){}

      void compute(PropagatorStub::WField const &)
      {}

   };

} 
#endif
