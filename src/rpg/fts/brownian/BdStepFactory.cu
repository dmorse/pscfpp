/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStepFactory.h"  
#include <rpg/fts/brownian/BdSimulator.h>

// Subclasses of BdStep 
#include <rpg/fts/brownian/ExplicitBdStep.h>
#include <rpg/fts/brownian/PredCorrBdStep.h>
#include <rpg/fts/brownian/LMBdStep.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   BdStepFactory<D>::BdStepFactory(Rp::BdSimulator<D, Rpg::Types<D> >& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of BdStep subclass className.
   */
   template <int D>
   Rp::BdStep<D, Rpg::Types<D> >* 
   BdStepFactory<D>::factory(const std::string &className) const
   {
      Rp::BdStep<D, Rpg::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      
      // Try to match classname
      if (className == "ExplicitBdStep" || className == "BdStep") {
         ptr = new Rp::ExplicitBdStep<D, Rpg::Types<D> >(*simulatorPtr_);
      } else
      if (className == "PredCorrBdStep") {
         ptr = new Rp::PredCorrBdStep<D, Rpg::Types<D> >(*simulatorPtr_);
      } else
      if (className == "LMBdStep") {
         ptr = new Rp::LMBdStep<D, Rpg::Types<D> >(*simulatorPtr_);
      }

      return ptr;
   }

   // Explicit instantiation definitions
   template class BdStepFactory<1>;
   template class BdStepFactory<2>;
   template class BdStepFactory<3>;

}
}
