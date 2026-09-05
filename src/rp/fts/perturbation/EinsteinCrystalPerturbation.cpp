/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>
#include <pscf/backend/cpp/CPT.h>

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EinsteinCrystalPerturbation<1,CPT>;
      template class EinsteinCrystalPerturbation<2,CPT>;
      template class EinsteinCrystalPerturbation<3,CPT>;
   }
}
