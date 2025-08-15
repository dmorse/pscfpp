#ifndef RPG_SIMULATOR_FACTORY_H
#define RPG_SIMULATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/fts/simulator/Simulator.h>

#include <string>

namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;

   /**
   * Factory for subclasses of Simulator.
   *
   * \ingroup Rpg_Fts_Module
   */
   template <int D>
   class SimulatorFactory : public Factory< Simulator<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System<D> object
      */
      SimulatorFactory(System<D>& system);

      /**
      * Method to create any Simulator supplied with PSCF.
      *
      * \param className name of the Simulator subclass
      * \return Simulator* pointer to new instance of className
      */
      Simulator<D>* factory(const std::string &className) const;

      using Factory< Simulator<D> >::trySubfactories;
      using Factory< Simulator<D> >::readObjectOptional;

   private:

      /// Pointer to the parent system.
      System<D>* systemPtr_;

   };

   #ifndef RPG_SIMULATOR_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class SimulatorFactory<1>;
   extern template class SimulatorFactory<2>;
   extern template class SimulatorFactory<3>;
   #endif

}
}
#endif
