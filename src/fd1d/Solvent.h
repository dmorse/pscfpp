#ifndef FD1D_SOLVENT_H
#define FD1D_SOLVENT_H

/*
* PFTS - Solvent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include <pfts/SolventTmpl.h>

namespace Pfts { 
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
