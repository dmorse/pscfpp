#ifndef PSPC_MC_MOVE_FACTORY_TPP
#define PSPC_MC_MOVE_FACTORY_TPP

#include "McMoveFactory.h"  

// Subclasses of McMove 
//#include "FourierMove.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   McMoveFactory<D>::McMoveFactory(McSimulator<D>& mcSimulator)
    : mcSimulatorPtr_(&mcSimulator)
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

      #if 0 
      // Try to match classname
      if (className == "FourierMove") {
         ptr = new FilmMcMove<D, AmMcMove<D> >(*mcSimulatorPtr_);
      }
      #endif

      return ptr;
   }

}
}
#endif
