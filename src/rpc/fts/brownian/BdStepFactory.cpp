/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStepFactory.h"
//#include <rpc/fts/brownian/BdSimulator.h>

// Subclasses of BdStep
#include "ExplicitBdStep.h"
#include "PredCorrBdStep.h"
#include "LMBdStep.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   BdStepFactory<D>::BdStepFactory(Rp::BdSimulator<D, Rpc::Types<D> >& simulator)
    : simulatorPtr_(&simulator)
   {}

   /*
   * Return a pointer to a instance of BdStep subclass className.
   */
   template <int D>
   Rp::BdStep<D, Rpc::Types<D> >* BdStepFactory<D>::factory(const std::string &className) const
   {
      Rp::BdStep<D, Rpc::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "ExplicitBdStep" || className == "BdStep") {
         ptr = new ExplicitBdStep<D>(*simulatorPtr_);
      } else
      if (className == "PredCorrBdStep") {
         ptr = new PredCorrBdStep<D>(*simulatorPtr_);
      } else
      if (className == "LMBdStep") {
         ptr = new LMBdStep<D>(*simulatorPtr_);
      }

      return ptr;
   }

   // Explicit instantiation definitions
   template class BdStepFactory<1>;
   template class BdStepFactory<2>;
   template class BdStepFactory<3>;

}
}
