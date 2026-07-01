/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationFactory.h"

// Subclasses of Perturbation
#include <rp/fts/perturbation/EinsteinCrystalPerturbation.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   PerturbationFactory<D,T>::PerturbationFactory(
		                Simulator<D,T>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /*
   * Return a pointer to a instance of Perturbation subclass className.
   */
   template <int D, class T>
   Perturbation<D,T>*
   PerturbationFactory<D,T>::factory(const std::string & className) const
   {
      Perturbation<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "EinsteinCrystal" ||
          className == "EinsteinCrystalPerturbation") {
         ptr = new EinsteinCrystalPerturbation<D,T>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
