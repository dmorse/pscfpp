#ifndef RP_MIXTURE_U_H
#define RP_MIXTURE_U_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureBase.h>   // base class template
#include <pscf/backend/cuda/CUT.h>    // template argument

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   template <int D, class T> class Mixture;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * Specializations of this template with D=1, 2, and 3 are derived
   * from specializations of base class template Rp::Mixture, and
   * inherit their public interface and almost all of their source
   * code from this base class.
   *
   * \see Rp::MixtureBase
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rp_Solver_Module
   */
   template <int D>
   class Mixture<D,CUT> 
     : public Rp::MixtureBase<D,CUT>
   {};

} // namespace Rp
} // namespace Pscf

#if 0
// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class MixtureTmpl< Rp::Polymer<1,CUT>, Rp::Solvent<1,CUT> >;
   extern template 
   class MixtureTmpl< Rp::Polymer<2,CUT>, Rp::Solvent<2,CUT> >;
   extern template 
   class MixtureTmpl< Rp::Polymer<3,CUT>, Rp::Solvent<3,CUT> >;
   namespace Rp {
      extern template class MixtureBase<1,CUT>;
      extern template class MixtureBase<2,CUT>;
      extern template class MixtureBase<3,CUT>;
   }
}
#endif

#endif
