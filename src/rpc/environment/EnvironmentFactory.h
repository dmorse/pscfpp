#ifndef RPC_ENVIRONMENT_FACTORY_H
#define RPC_ENVIRONMENT_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/environment/Environment.h>
#include <rpc/system/Types.h>
#include <util/param/Factory.h>  

#include <string>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rpc {


   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Factory for subclasses of Environment.
   *
   * \ingroup Rpc_Environment_Module
   */

   template <int D>
   class EnvironmentFactory : public Factory<Environment> 
   {

   public:

      /**
      * Constructor.
      */
      EnvironmentFactory(Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Method to create any Environment supplied with PSCF.
      *
      * \param className name of the Environment subclass
      * \return Environment* pointer to new instance of className
      */
      Environment* factory(const std::string &className) const;

      using Factory<Environment>::trySubfactories;
      using Factory<Environment>::readObjectOptional;

   private:

      /// Pointer to the parent system.
      Rp::System<D, Rpc::Types<D> >* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class EnvironmentFactory<1>;
   extern template class EnvironmentFactory<2>;
   extern template class EnvironmentFactory<3>;

}
}
#endif
