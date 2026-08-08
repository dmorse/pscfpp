/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimulatorFactory.h"
#include <rp/system/System.h>

#include <rp/fts/brownian/BdSimulator.h>
//#include <rp/fts/montecarlo/McSimulator.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   SimulatorFactory<D,T>::SimulatorFactory(System<D,T>& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Simulator subclass className.
   */
   template <int D, class T>
   Simulator<D,T>* SimulatorFactory<D,T>::factory(const std::string &className)
   const
   {
      Simulator<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "BdSimulator" || className == "Simulator") {
         ptr = new BdSimulator<D,T>(*systemPtr_);
      } else
      if (className == "McSimulator") {
         ptr = new McSimulator<D,T>(*systemPtr_);
      } 

      return ptr;
   }

}
}
