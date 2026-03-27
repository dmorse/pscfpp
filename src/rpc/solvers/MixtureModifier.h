#ifndef RPC_MIXTURE_MODIFIER_H
#define RPC_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.h>  // base class template
#include "Mixture.h"                     // base class argument

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Parameter modifier for an associated Mixture.
   *
   * Each specialization of this template is simply a named instantation 
   * of the base class template Rp::MixtureModifier and has the same public
   * interface as this base class. 
   *
   * \see Rp::MixtureModifier
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class MixtureModifier : public Rp::MixtureModifier< Mixture<D> >
   {};

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class MixtureModifier< Rpc::Mixture<1> >;
      extern template class MixtureModifier< Rpc::Mixture<2> >;
      extern template class MixtureModifier< Rpc::Mixture<3> >;
   }
   namespace Rpc {
      extern template class MixtureModifier<1>;
      extern template class MixtureModifier<2>;
      extern template class MixtureModifier<3>;
   }        
} 
#endif
