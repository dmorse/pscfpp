/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1,CUT>;
      template class ScftThermo<2,CUT>;
      template class ScftThermo<3,CUT>;
   }
}
