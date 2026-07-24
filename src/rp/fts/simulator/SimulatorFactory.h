#ifndef RP_SIMULATOR_FACTORY_H
#define RP_SIMULATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>          // base class template
#include <rp/fts/simulator/Simulator.h>  // base class argument
#include <pscf/backends/TmplDeclare.h>   // declaration macros

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

   /**
   * Factory for subclasses of Simulator.
   *
   * \ingroup Rp_Fts_Simulator_Module
   */
   template <int D, class T>
   class SimulatorFactory : public Factory< Simulator<D,T> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System<D,T> object
      */
      SimulatorFactory(System<D,T>& system);

      /**
      * Method to create any Simulator supplied with PSCF.
      *
      * \param className name of the Simulator subclass
      * \return Simulator* pointer to new instance of className
      */
      Simulator<D,T>* factory(const std::string &className) const;

      using Factory< Simulator<D,T> >::trySubfactories;
      using Factory< Simulator<D,T> >::readObjectOptional;

   private:

      /// Pointer to the parent system.
      System<D,T>* systemPtr_;

   };

   // Explicit instantation declarations
   PSCF_TMPL_DECLARE(SimulatorFactory)

}
}
#endif
