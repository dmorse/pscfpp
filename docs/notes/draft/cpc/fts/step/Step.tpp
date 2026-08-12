#ifndef CPC_BD_STEP_TPP
#define CPC_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Step.h"

#include <cpc/fts/Simulator.h>
#include <cpc/system/System.h>

namespace Pscf {
namespace Cpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Step<D>::Step(Simulator<D>& simulator)
    : simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      randomPtr_(&(simulator.random()))
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   Step<D>::~Step()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void Step<D>::readParameters(std::istream &in)
   {}

   /*
   * Setup at beginning of loop.
   */
   template <int D>
   void Step<D>::setup()
   {}

   template <int D>
   void Step<D>::output()
   {}

   template<int D>
   void Step<D>::outputTimers(std::ostream& out)
   {}

   template<int D>
   void Step<D>::clearTimers()
   {}

}
}
#endif
