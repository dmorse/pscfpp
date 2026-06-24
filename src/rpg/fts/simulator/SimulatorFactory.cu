/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimulatorFactory.h"  
#include <rpg/system/System.h>

// Subclasses of Simulator 
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/brownian/BdSimulator.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   SimulatorFactory<D>::SimulatorFactory(Rp::System<D, Rpg::Types<D> >& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Simulator subclass className.
   */
   template <int D>
   Rp::Simulator<D, Rpg::Types<D> >* SimulatorFactory<D>::factory(const std::string &className) 
   const
   {
      Rp::Simulator<D, Rpg::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
      
      // Try to match classname
      if (className == "McSimulator" || className == "Simulator") {
         ptr = new Rp::McSimulator<D, Rpg::Types<D> >(*systemPtr_);
      } else
      if (className == "BdSimulator") {
         ptr = new BdSimulator<D>(*systemPtr_);
      } 

      return ptr;
   }


   // Explicit instantiation definitions
   template class SimulatorFactory<1>;
   template class SimulatorFactory<2>;
   template class SimulatorFactory<3>;

}
}
