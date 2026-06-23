/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimulatorFactory.h"
//#include <rpc/system/System.h>

// Subclasses of Simulator
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/fts/brownian/BdSimulator.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   SimulatorFactory<D>::SimulatorFactory(Rp::System<D, Rpc::Types<D> >& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Simulator subclass className.
   */
   template <int D>
   Simulator<D>* SimulatorFactory<D>::factory(const std::string &className)
   const
   {
      Simulator<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "McSimulator" || className == "Simulator") {
         ptr = new McSimulator<D>(*systemPtr_);
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
