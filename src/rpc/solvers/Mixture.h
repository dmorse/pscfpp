#ifndef RPC_MIXTURE_H
#define RPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Mixture.h>   // base class template
#include <rpc/system/Types.h>     // base class template argument

namespace Pscf {
namespace Rpc {

   // Forward declarations
   template <int D> class Polymer;
   template <int D> class Solvent;

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of the base class template Rp::Mixture,
   * and inherit their public interface and almost all of their source
   * code from this base class.
   *
   * \see Rp::Mixture
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Mixture : public Rp::Mixture< D, Types<D> >
   {
   public:

      /// Direct base class
      using RpMixtureT = typename Rp::Mixture< D, Types<D> >;

      // Inherited names
      using typename RpMixtureT::MixtureBaseT;
      using typename RpMixtureT::MixtureTmplT;
      using typename RpMixtureT::FieldT;
      using MixtureTmplT::polymer;

   private:

      /**
      * Allocate memory for all blocks
      */
      void allocateBlocks() override;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   extern template class MixtureTmpl< Rpc::Polymer<1>, Rpc::Solvent<1> >;
   extern template class MixtureTmpl< Rpc::Polymer<2>, Rpc::Solvent<2> >;
   extern template class MixtureTmpl< Rpc::Polymer<3>, Rpc::Solvent<3> >;
   namespace Rp {
      extern template class Mixture<1, Rpc::Types<1> >;
      extern template class Mixture<2, Rpc::Types<2> >;
      extern template class Mixture<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Mixture<1>;
      extern template class Mixture<2>;
      extern template class Mixture<3>;
   }
}
#endif
