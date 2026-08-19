/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RampFactory.h"  

// Subclasses of Ramp 
#include "LinearRamp.h"

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   RampFactory<D,T>::RampFactory(Simulator<D,T>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of Ramp subclass className.
   */
   template <int D, class T>
   Ramp<D, T >* RampFactory<D,T>::factory(const std::string & className) const
   {
      Ramp<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
       
      // Try to match classname
      if (className == "LinearRamp") {
         ptr = new LinearRamp<D,T>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
