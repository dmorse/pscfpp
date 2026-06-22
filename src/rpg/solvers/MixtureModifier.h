#ifndef RPG_MIXTURE_MODIFIER_H
#define RPG_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.h>  // base class template
#include <rpg/solvers/Mixture.h>         // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Prdc;

   /**
   * Parameter modifier for an associated Mixture.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from corresponding specializations of the base class template 
   * Rp::MixtureModifier, and inherit their public interface and almost 
   * all of their source code from this base class.  
   *
   * \see Rp::MixtureModifier
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class MixtureModifier : public Rp::MixtureModifier< Rp::Mixture<D, Rpg::Types<D> > >
   { 
   public:

      MixtureModifier() = default;

      virtual ~MixtureModifier() = default;
   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class MixtureModifier< Rp::Mixture<1, Rpg::Types<1> > >;
      extern template class MixtureModifier< Rp::Mixture<2, Rpg::Types<2> > >;
      extern template class MixtureModifier< Rp::Mixture<3, Rpg::Types<3> > >;
   } 
   namespace Rpg {
      extern template class MixtureModifier<1>;
      extern template class MixtureModifier<2>;
      extern template class MixtureModifier<3>;
   } 
}
#endif
