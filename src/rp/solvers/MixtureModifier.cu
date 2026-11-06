/*
* PSCF - MixtureModifier Self-Consistent Field Theory
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MixtureModifier<1,CUT>;
      template class MixtureModifier<2,CUT>;
      template class MixtureModifier<3,CUT>;
   }
}
