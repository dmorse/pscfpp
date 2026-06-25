/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationFactory.h"  
//#include <rpc/fts/simulator/Simulator.h>

// Subclasses of Perturbation 
#include "EinsteinCrystalPerturbation.h"


namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   PerturbationFactory<D>::PerturbationFactory(Rp::Simulator<D, Rpc::Types<D> >& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of Perturbation subclass className.
   */
   template <int D>
   Rp::Perturbation<D, Rpc::Types<D> >* 
   PerturbationFactory<D>::factory(const std::string & className) const
   {
      Rp::Perturbation<D, Rpc::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
       
      // Try to match classname
      if (className == "EinsteinCrystal" || 
          className == "EinsteinCrystalPerturbation") {
         ptr = new EinsteinCrystalPerturbation<D>(*simulatorPtr_);
      } 

      return ptr;
   }

   // Explicit instantiation definitions
   template class PerturbationFactory<1>;
   template class PerturbationFactory<2>;
   template class PerturbationFactory<3>;

}
}
