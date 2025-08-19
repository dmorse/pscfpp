#ifndef RPC_MIXTURE_MODIFIER_H
#define RPC_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/solvers/MixtureModifierPrdc.h>  // base class template
#include "Mixture.h"                           // base class parameter

namespace Pscf {
namespace Rpc {

   /**
   * Parameter modifier for an associated Mixture.
   *
   * A MixtureModifier<D> is derived from the template specialization 
   * Prdc::MixtureModifierPrdc< Rpc::Mixture<D> > and has the same public
   * interface as this base class. See documentation of base class.
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class MixtureModifier : public MixtureModifierPrdc< Mixture<D> >
   {

   public:
 
      /// Direct (parent) base class.
      using Base = typename Prdc::MixtureModifierPrdc< Mixture<D> >;

      // Inherited public member functions

      using Base::associate;
      using Base::setKuhn;
      using Base::setPhiPolymer;
      using Base::setMuPolymer;
      using Base::setBlockLength;
      using Base::setPhiSolvent;
      using Base::setMuSolvent;
      using Base::setSolventSize;
      using Base::setVMonomer;

   };

   // Suppress implicit instantiation
   extern template class MixtureModifier<1>;
   extern template class MixtureModifier<2>;
   extern template class MixtureModifier<3>;

} // namespace Rpc
namespace Prdc {
   // Suppress implicit instantiation of base class 
   extern template class MixtureModifierPrdc< Rpc::Mixture<1> >;
   extern template class MixtureModifierPrdc< Rpc::Mixture<2> >;
   extern template class MixtureModifierPrdc< Rpc::Mixture<3> >;
} // namespace Prdc
} // namespace Pscf
#endif
