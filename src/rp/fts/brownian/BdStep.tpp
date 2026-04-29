#ifndef RP_BD_STEP_TPP
#define RP_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStep.h"
#include <util/random/Random.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   BdStep<D,T>::BdStep(typename T::BdSimulator& simulator)
    : simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      randomPtr_(&(simulator.random())),
      vecRandomPtr_(&(simulator.vecRandom()))
   {}

   /*
   * readParameters, empty default implementation.
   */
   template <int D, class T>
   void BdStep<D,T>::readParameters(std::istream &in)
   {}

   /*
   * Setup at beginning of loop.
   */
   template <int D, class T>
   void BdStep<D,T>::setup()
   {}

   template <int D, class T>
   void BdStep<D,T>::output()
   {}

   template <int D, class T>
   void BdStep<D,T>::outputTimers(std::ostream& out)
   {}

   template<int D, class T>
   void BdStep<D,T>::clearTimers()
   {}

}
}
#endif
