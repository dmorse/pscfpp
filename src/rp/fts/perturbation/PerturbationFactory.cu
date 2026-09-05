/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>
#include <rp/fts/perturbation/PerturbationFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationFactory<1,CUT>;
      template class PerturbationFactory<2,CUT>;
      template class PerturbationFactory<3,CUT>;
   }
}
