#ifndef RPG_RAMP_TPP
#define RPG_RAMP_TPP

#include "Ramp.h"
#include <rpg/fts/Simulator.h>
#include <prdc/cuda/RField.h>

#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   Ramp<D>::Ramp(Simulator<D>& simulator)
    : ParamComposite(),
      simulatorPtr_(&simulator)
   {}
   
   /* 
   * Destructor.
   */
   template <int D>
   Ramp<D>::~Ramp()
   {}

   /* 
   * Setup before simulation - sets the nStep member variable. 
   */
   template <int D>
   void Ramp<D>::setup(int nStep)
   {  nStep_ = nStep; }

}
}
#endif 
