#ifndef PSCF_SOLVENT_STUB_H
#define PSCF_SOLVENT_STUB_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>
#include <util/param/ParamComposite.h>

namespace Pscf
{ 

   class SolventStub : public SolventSpecies<double>
   {

   public:

      SolventStub()
      {  setClassName("Solvent"); }

      void compute(PropagatorStub::FieldT const &)
      {}

   };

} 
#endif
