#ifndef PSCF_SOLVENT_STUB_H
#define PSCF_SOLVENT_STUB_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>
#include <util/param/ParamComposite.h>

namespace Pscf
{ 

   class SolventStub : public SolventSpecies
   {

   public:

      SolventStub()
      {  setClassName("Solvent"); }

      void compute(PropagatorStub::WField const &)
      {}

   };

} 
#endif
