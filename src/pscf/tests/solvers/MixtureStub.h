#ifndef PSCF_MIXTURE_STUB_H
#define PSCF_MIXTURE_STUB_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerStub.h"
#include "SolventStub.h"
#include <pscf/solvers/MixtureTmpl.h>

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
