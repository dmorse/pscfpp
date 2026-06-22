#ifndef RPG_SYSTEM_H
#define RPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <rp/system/System.h>       // base class template
#include <rpg/system/Types.h>       // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Main class, representing a complete physical system.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specialization of the base class template Rp::System,
   * and inherit their entire public interface and almost all of their
   * source code from this base class.
   *
   * \see Rp::System
   * \ingroup Rpg_System_Module
   */
   template <int D>
   class System : public Rp::System< D, Types<D> >
   {
   public:

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      virtual ~System() = default;

      // Prohibit copying and assignment.
      System(System<D> const &) = delete;
      System<D>& operator = (System<D> const &) = delete;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::System<1, Rpg::Types<1> >;
      extern template class Rp::System<2, Rpg::Types<1> >;
      extern template class Rp::System<3, Rpg::Types<1> >;
   }
   namespace Rpg {
      extern template class System<1>;
      extern template class System<2>;
      extern template class System<3>;
   }
}
#endif
