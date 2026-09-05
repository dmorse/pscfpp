/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>
#include <pscf/backend/cuda/CUT.h>

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EinsteinCrystalPerturbation<1,CUT>;
      template class EinsteinCrystalPerturbation<2,CUT>;
      template class EinsteinCrystalPerturbation<3,CUT>;
   }
}
