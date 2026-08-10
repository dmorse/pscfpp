/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
#include <rp/fts/perturbation/PerturbationFactory.h>

// Subclasses of Perturbation
#include <rp/fts/perturbation/EinsteinCrystalPerturbation.h>
#endif

#include <rp/fts/perturbation/PerturbationFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationFactory<1,CUT>;
      template class PerturbationFactory<2,CUT>;
      template class PerturbationFactory<3,CUT>;
   }
}
