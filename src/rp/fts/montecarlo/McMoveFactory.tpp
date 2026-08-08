/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMoveFactory.h"
//#include <rp/fts/montecarlo/McSimulator.h>

// Subclasses of McMove
#include <rp/fts/montecarlo/RealMove.h>
//#include <rp/fts/montecarlo/ForceBiasMove.h>
#include <rp/fts/montecarlo/BdMove.h>
//#include <rp/fts/montecarlo/ShiftMove.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   McMoveFactory<D,T>::McMoveFactory(McSimulator<D,T>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /*
   * Return a pointer to a instance of McMove subclass className.
   */
   template <int D, class T>
   McMove<D,T>* 
   McMoveFactory<D,T>::factory(const std::string &className) const
   {
      McMove<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "RealMove") {
         ptr = new RealMove<D,T>(*simulatorPtr_);
      } else if (className == "ForceBiasMove") {
         ptr = new ForceBiasMove<D,T>(*simulatorPtr_);
      } else if (className == "BdMove") {
         ptr = new BdMove<D,T>(*simulatorPtr_);
      } else if (className == "ShiftMove") {
         ptr = new ShiftMove<D,T>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
