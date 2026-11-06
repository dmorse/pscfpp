/*
* PSCF - MixtureModifier Self-Consistent Field Theory
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions (CPU)
namespace Pscf {
   namespace Rp { 
      template class MixtureModifier<1,CPT>;
      template class MixtureModifier<2,CPT>;
      template class MixtureModifier<3,CPT>;
   }
}
