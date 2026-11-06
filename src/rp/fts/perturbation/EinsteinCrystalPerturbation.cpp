/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EinsteinCrystalPerturbation<1,CPT>;
      template class EinsteinCrystalPerturbation<2,CPT>;
      template class EinsteinCrystalPerturbation<3,CPT>;
   }
}
