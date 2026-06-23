#ifndef RPG_SYSTEM_CONST_REF_H
#define RPG_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>   // base class template
#include <rpg/system/Types.h>           // base class parameter

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Const access to a System.
   *
   * Specializations of this template with D=1, 2, and 3 are derived
   * from corresponding specializations of the base class template
   * Rp::SystemConstRef, and inherit their public interface and almost
   * all of their source code from this base class.
   *
   * \see Rp::SystemConstRef
   * \ingroup Rpg_System_Module
   */
   template <int D>
   class SystemConstRef : public Rp::SystemConstRef<D, Types<D> >
   {
   public:

      /**
      * Default constructor.
      */
      SystemConstRef()
       : Rp::SystemConstRef< D, Types<D> > ()
      {};

      /**
      * Constructor.
      *
      * \param system  Rp::System<D, Rpg::Types<D> > object to which this refers.
      */
      SystemConstRef(Rp::System<D, Rpg::Types<D> > const & system)
       : Rp::SystemConstRef<D, Types<D> > (system)
      {};

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef<1, Rpg::Types<1> >;
      extern template class SystemConstRef<2, Rpg::Types<2> >;
      extern template class SystemConstRef<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class SystemConstRef<1>;
      extern template class SystemConstRef<2>;
      extern template class SystemConstRef<3>;
   }
}
#endif
