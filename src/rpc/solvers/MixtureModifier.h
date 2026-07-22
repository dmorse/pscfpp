#ifndef RPC_MIXTURE_MODIFIER_H
#define RPC_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.h>  // base class template
#include <pscf/backends/CPT.h>            // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Pscf::Prdc;
      extern template class MixtureModifier<1,CPT>;
      extern template class MixtureModifier<2,CPT>;
      extern template class MixtureModifier<3,CPT>;
   }
} 
#endif
