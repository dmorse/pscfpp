#ifndef PSCF_SOLVENT_STUB_H
#define PSCF_SOLVENT_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/SolventTmpl.h>
#include <util/containers/DArray.h>

namespace Pscf
{ 

   class SolventStub : public SolventTmpl< DArray<double> >
   {

   public:

      SolventStub(){}
   
   };

} 
#endif
