#ifndef CHEM_SOLVENT_STUB_H
#define CHEM_SOLVENT_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <chem/SolventTmpl.h>
#include <util/containers/DArray.h>

namespace Chem{ 

   class SolventStub : public SolventTmpl< DArray<double> >
   {

   public:

      SolventStub(){}
   
   };

} 
#endif
