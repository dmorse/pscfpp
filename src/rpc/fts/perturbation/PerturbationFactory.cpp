/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/perturbation/PerturbationFactory.h>

// Subclasses of Perturbation
#include <rpc/fts/perturbation/EinsteinCrystalPerturbation.h>

#include <rp/fts/perturbation/PerturbationFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationFactory<1, CppTp<1> >;
      template class PerturbationFactory<2, CppTp<2> >;
      template class PerturbationFactory<3, CppTp<3> >;
   }
}
