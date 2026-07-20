#ifndef RPC_PERTURBATION_FACTORY_H
#define RPC_PERTURBATION_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/perturbation/PerturbationFactory.h>
#include <pscf/cpu/CppTp.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PerturbationFactory<1, CppTp<1> >;
      extern template class PerturbationFactory<2, CppTp<2> >;
      extern template class PerturbationFactory<3, CppTp<3> >;
   }
}
#endif
