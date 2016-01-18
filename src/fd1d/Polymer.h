#ifndef FD1D_POLYMER_H
#define FD1D_POLYMER_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/PolymerTmpl.h>

namespace Pscf { 
namespace Fd1d
{ 

   class Polymer : public PolymerTmpl<Block>
   {

   public:

      Polymer();

      ~Polymer();

   };

} 
}
#endif
