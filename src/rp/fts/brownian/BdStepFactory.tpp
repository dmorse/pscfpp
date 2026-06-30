/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStepFactory.h"
//#include <rp/fts/brownian/BdSimulator.h>

// Subclasses of BdStep
#include <rp/fts/brownian/ExplicitBdStep.h>
#include <rp/fts/brownian/PredCorrBdStep.h>
#include <rp/fts/brownian/LMBdStep.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc; 

   /*
   * Constructor
   */
   template <int D, class T>
   BdStepFactory<D,T>::BdStepFactory(BdSimulator<D,T>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /*
   * Return a pointer to a instance of BdStep subclass className.
   */
   template <int D, class T>
   BdStep<D,T>* BdStepFactory<D,T>::factory(const std::string &className) const
   {
      BdStep<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "ExplicitBdStep" || className == "BdStep") {
         ptr = new ExplicitBdStep<D,T>(*simulatorPtr_);
      } else
      if (className == "PredCorrBdStep") {
         ptr = new PredCorrBdStep<D,T>(*simulatorPtr_);
      } else
      if (className == "LMBdStep") {
         ptr = new LMBdStep<D,T>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
