#ifndef RPC_PERTURBATION_FACTORY_TPP
#define RPC_PERTURBATION_FACTORY_TPP

#include "PerturbationFactory.h"  

// Subclasses of Perturbation 
//#include "EinsteinCrystalPerturbation.h"

#include <rpc/simulate/Simulator.h>

namespace Pscf {
namespace Rpc {

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
      //if (className == "EinsteinCrystal") {
      //   ptr = new EinsteinCrystal<D>(*simulatorPtr_);
      //} 

      return ptr;
   }

}
}
#endif
