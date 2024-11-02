#ifndef RPG_PERTURBATION_FACTORY_TPP
#define RPG_PERTURBATION_FACTORY_TPP

#include "PerturbationFactory.h"  

// Subclasses of Perturbation 
#include "EinsteinCrystalPerturbation.h"

#include <rpg/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   PerturbationFactory<D>::PerturbationFactory(Simulator<D>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of Perturbation subclass className.
   */
   template <int D>
   Perturbation<D>* 
   PerturbationFactory<D>::factory(const std::string & className) const
   {
      Perturbation<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
       
      // Try to match classname
      if (className == "EinsteinCrystal" || 
          className == "EinsteinCrystalPerturbation") {
         ptr = new EinsteinCrystalPerturbation<D>(*simulatorPtr_);
      } 

      return ptr;
   }

}
}
#endif
