#ifndef FD1D_SOLVENT_H
#define FD1D_SOLVENT_H

/*
* PFTS - Solvent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include <pscf/SolventTmpl.h>

namespace Pscf { 
namespace Fd1d
{ 

   class Solvent : public SolventTmpl<Propagator>
   {
   public:

      /**
      * Constructor.
      */
      Solvent();

      /**
      * Destructor.
      */
      ~Solvent();

   };

}
} 
#endif
