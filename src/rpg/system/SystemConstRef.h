#ifndef RPG_SYSTEM_CONST_REF_H
#define RPG_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/system/SystemConstRefReal.h>   // base class template
#include <rpg/System.h>                       // template parameter

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Const access to a System<D>.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class SystemConstRef : public SystemConstRefReal< System<D> >
   {
   public:

      /// Alias for base class
      using Base = SystemConstRefReal< System<D> >;

      /**
      * Default constructor.
      */
      SystemConstRef()
       : Base()
      {};

      /**
      * Constructor.
      */
      SystemConstRef(System<D> const & system)
       : Base(system)
      {};

   };

   // Suppress implicit instantiation
   extern template class SystemConstRef<1>;
   extern template class SystemConstRef<2>;
   extern template class SystemConstRef<3>;

}

namespace Prdc {
   // Suppress implicit instantiation of base class
   extern template class SystemConstRefReal< Rpg::System<1> >;
   extern template class SystemConstRefReal< Rpg::System<2> >;
   extern template class SystemConstRefReal< Rpg::System<3> >;
}

} // namespace Pscf
#endif
