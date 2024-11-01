#ifndef RPG_MC_MOVE_FACTORY_TPP
#define RPG_MC_MOVE_FACTORY_TPP

#include "McMoveFactory.h"  
#include <rpg/fts/montecarlo/McSimulator.h>

// Subclasses of McMove 
#include "RealMove.h"
#include "FourierMove.h"
#include "ForceBiasMove.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   McMoveFactory<D>::McMoveFactory(McSimulator<D>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of McMove subclass className.
   */
   template <int D>
   McMove<D>* McMoveFactory<D>::factory(const std::string &className) const
   {
      McMove<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      
      // Try to match classname
      if (className == "RealMove") {
         ptr = new RealMove<D>(*simulatorPtr_);
      } else if (className == "FourierMove") {
         ptr = new FourierMove<D>(*simulatorPtr_);
      } else if (className == "ForceBiasMove") {
         ptr = new ForceBiasMove<D>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
#endif
