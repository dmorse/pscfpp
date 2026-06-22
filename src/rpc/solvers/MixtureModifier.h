#ifndef RPC_MIXTURE_MODIFIER_H
#define RPC_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.h>  // base class template
#include <rpc/solvers/Mixture.h>         // base class argument

#if 0
namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Parameter modifier for an associated Mixture.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of base class template Rp::MixtureModifier, 
   * and inherit their public interface and all of their source code
   * from this base class.  
   *
   * \see Rp::MixtureModifier
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class MixtureModifier : public Rp::MixtureModifier< Rp::Mixture<D, Rpc::Types<D> > >
   {};

} // namespace Rpc
} // namespace Pscf
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class MixtureModifier<1, Rpc::Types<1> >;
      extern template class MixtureModifier<2, Rpc::Types<2> >;
      extern template class MixtureModifier<3, Rpc::Types<3> >;
   }
} 
#endif
