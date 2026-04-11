#ifndef RPG_SYSTEM_CONST_REF_H
#define RPG_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>  // base class template
#include <rpg/system/System.h>         // template parameter argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * Const access to a System<D>.
   *
   * Specializations of this template with D=1, 2, and 3 are derived
   * from specializations of base class template Rp::SystemConstRef,
   * and inherit their public interface and almost all of their source
   * code from this base class.
   *
   * \see Rp::SystemConstRef
   * \ingroup Rpg_System_Module
   */
   template <int D>
   class SystemConstRef : public Rp::SystemConstRef< System<D> >
   {
   public:

      /// Alias for base class.
      using RpSystemConstRef = Rp::SystemConstRef< System<D> >;

      /**
      * Default constructor.
      */
      SystemConstRef()
       : RpSystemConstRef()
      {};

      /**
      * Constructor (creates association with parent System).
      *
      * \param system  parent system
      */
      SystemConstRef(System<D> const & system)
       : RpSystemConstRef(system)
      {};

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef< Rpg::System<1> >;
      extern template class SystemConstRef< Rpg::System<2> >;
      extern template class SystemConstRef< Rpg::System<3> >;
   }
   namespace Rpg {
      extern template class SystemConstRef<1>;
      extern template class SystemConstRef<2>;
      extern template class SystemConstRef<3>;
   }
}
#endif
