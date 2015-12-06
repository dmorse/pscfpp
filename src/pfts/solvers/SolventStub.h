#ifndef PFTS_SOLVENT_STUB_H
#define PFTS_SOLVENT_STUB_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/solvers/SolventTmpl.h>
#include <util/containers/DArray.h>

namespace Pfts{ 

   class SolventStub : public SolventTmpl< DArray<double> >
   {

   public:

      SolventStub();
   
   };

} 
#endif
