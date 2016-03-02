#ifndef PSCF_MIXTURE_STUB_H
#define PSCF_MIXTURE_STUB_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerStub.h"
#include "SolventStub.h"
#include <pscf/MixtureTmpl.h>

namespace Pscf
{ 

   class MixtureStub 
    : public MixtureTmpl<PolymerStub, SolventStub>
   {

   public:

      MixtureStub()
      {  setClassName("Mixture"); }
   
   };

} 
#endif
