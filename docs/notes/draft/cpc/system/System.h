#ifndef CPC_SYSTEM_H
#define CPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <cp/system/System.h>            // base class template
#include <cpc/system/Types.h>     // base class template argument
#include <cpc/field/WFields.h>    // member of base class
#include <cpc/field/CFields.h>    // member of base class

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Main class for CL-FTS, representing a complete physical system.
   *
   * This class is derived from a specialization of the class template
   * Cp::System, and has the same public interface as this base 
   * class.  See the documentation of this base class template for 
   * details.
   *
   * \ingroup Cpc_System_Module
   */
   template <int D>
   class System : public Cp::System< D, Types<D> >
   {
   public:

      /// Inherit default constructor.
      using Cp::System< D, Types<D> >::System;

      /// Copy constructor.
      System(System<D> const &) = delete;

      /// Destructor.
      ~System() = default;

      /// Assignment
      System<D>& operator = (System<D> const &) = delete;

   };

} // namespace Cpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Cp {
      extern template class System<1, Cpc::Types<1> >;
      extern template class System<2, Cpc::Types<1> >;
      extern template class System<3, Cpc::Types<1> >;
   }
   namespace Cpc {
      extern template class System<1>;
      extern template class System<2>;
      extern template class System<3>;
   }
}
#endif
