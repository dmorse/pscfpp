#ifndef RPG_BASIS_FIELD_STATE_H
#define RPG_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/BasisFieldState.h>
#include <rpg/system/Types.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of base class template Rp::BasisFieldState, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::BasisFieldState
   * \ingroup Rpg_Scft_Sweep_Module
   */
   template <int D>
   class BasisFieldState
    : public Rp::BasisFieldState<D, Types<D> >
   {
   public:

      /**
      * Default constructor.
      */
      BasisFieldState();

      /**
      * Constructor, create association with a parent system.
      *
      * \param system associated parent system
      */
      BasisFieldState(Rp::System<D, Rpg::Types<D> >& system);

      virtual ~BasisFieldState() = default;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BasisFieldState< 1, Rpg::Types<1> >;
      extern template class BasisFieldState< 2, Rpg::Types<2> >;
      extern template class BasisFieldState< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BasisFieldState<1>;
      extern template class BasisFieldState<2>;
      extern template class BasisFieldState<3>;
   }
}
#endif
