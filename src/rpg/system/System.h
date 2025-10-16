#ifndef RPG_SYSTEM_H
#define RPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <prdc/system/SystemTmpl.h>    // base class template
#include <rpg/system/Types.h>          // base class template param
#include <rpg/field/WFields.h>         // member
#include <rpg/field/CFields.h>         // member
#include <rpg/field/Mask.h>            // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Main class, representing a complete physical system.
   *
   * This class is derived from a partial specialization of the base class
   * template Prdc::SystemTmpl, and has the same public interface as the
   * base class. See the documentation of this template for details.
   *
   * \ingroup Rpg_System_Module
   */
   template <int D>
   class System : public SystemTmpl< D, Types<D> >
   {
   public:

      /**
      * Constructor.
      */
      System();

      System(System<D> const &) = delete;
      System<D>& operator = (System<D> const &) = delete;

   protected:

      /**
      * Explicitly set maximum number of threads per block.
      *
      * This function is called in the setOptions function if the
      * -t command line option is present. 
      *
      * \param nThread  thread count (maximum per block). 
      */
      virtual void setThreadCount(int nThread) override;

   };

   // Explicit instantiation declarations
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;

} // namespace Rpg
namespace Prdc {

   // Explicit instantiation declarations of base class
   extern template class SystemTmpl<1, Rpg::Types<1> >;
   extern template class SystemTmpl<2, Rpg::Types<1> >;
   extern template class SystemTmpl<3, Rpg::Types<1> >;

} // namespace Prdc 
} // namespace Pscf
#endif
