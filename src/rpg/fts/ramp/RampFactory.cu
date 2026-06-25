/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RampFactory.h"  
#include <rpg/fts/simulator/Simulator.h>

// Subclasses of Ramp 
#include "LinearRamp.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   RampFactory<D>::RampFactory(Rp::Simulator<D, Rpg::Types<D> >& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of Ramp subclass className.
   */
   template <int D>
   Rp::Ramp<D, Rpg::Types<D> >* RampFactory<D>::factory(const std::string & className) const
   {
      Rp::Ramp<D, Rpg::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
       
      // Try to match classname
      if (className == "LinearRamp") {
         ptr = new Rp::LinearRamp<D, Rpg::Types<D> >(*simulatorPtr_);
      }

      return ptr;
   }


   // Explicit instantiation definitions
   template class RampFactory<1>;
   template class RampFactory<2>;
   template class RampFactory<3>;

}
}
