#ifndef RP_ENVIRONMENT_FACTORY_H
#define RP_ENVIRONMENT_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>            // base class template
#include <prdc/environment/Environment.h>  // base class argument
#include <pscf/backends/TmplDeclare.h>     // preprocessor macros

#include <string>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
   }
}

namespace Pscf {
namespace Rp {


   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Factory for subclasses of Environment.
   *
   * \ingroup Rp_Environment_Module
   */

   template <int D, class T>
   class EnvironmentFactory : public Factory<Environment> 
   {

   public:

      /**
      * Constructor.
      */
      EnvironmentFactory(System<D,T>& system);

      /**
      * Method to create any Environment supplied with PSCF.
      *
      * \param className name of the Environment subclass
      * \return Environment* pointer to new instance of className
      */
      Environment* factory(const std::string &className) const;

      //using Factory<Environment>::trySubfactories;
      //using Factory<Environment>::readObjectOptional;

   private:

      /// Pointer to the parent system.
      System<D,T>* sysPtr_;

   };

   // Explicit instantation declarations
   PSCF_TMPL_DECLARE(EnvironmentFactory)

}
}
#endif
