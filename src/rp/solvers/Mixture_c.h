#ifndef RP_MIXTURE_C_H
#define RP_MIXTURE_C_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureBase.h>   // base class template
#include <pscf/backend/cpp/CPT.h>     // template argument

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   template <int D, class T> class Mixture;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * This is the partial specialization for use with a CPU backend.
   * The public interface is identical to that of the base class.
   *
   * \see Rp::MixtureBase
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rp_Solver_Module
   */
   template <int D>
   class Mixture<D,CPT> 
    : public Rp::MixtureBase<D,CPT>
   {};

} // namespace Rp
} // namespace Pscf

#if 0
// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class MixtureTmpl< Rp::Polymer<1,CPT>, Rp::Solvent<1,CPT> >;
   extern template 
   class MixtureTmpl< Rp::Polymer<2,CPT>, Rp::Solvent<2,CPT> >;
   extern template 
   class MixtureTmpl< Rp::Polymer<3,CPT>, Rp::Solvent<3,CPT> >;
   namespace Rp {
      extern template class MixtureBase<1,CPT>;
      extern template class MixtureBase<2,CPT>;
      extern template class MixtureBase<3,CPT>;
   }
}
#endif

#endif
