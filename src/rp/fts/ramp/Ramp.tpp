#ifndef RP_RAMP_TPP
#define RP_RAMP_TPP

#include "Ramp.h"

#include <rp/fts/simulator/Simulator.h>
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D, class T>
   Ramp<D,T>::Ramp(Simulator<D,T>& simulator)
    : ParamComposite(),
      simulatorPtr_(&simulator)
   {}
   
   /* 
   * Setup before simulation - sets the nStep member variable. 
   */
   template <int D, class T>
   void Ramp<D,T>::setup(int nStep)
   {  nStep_ = nStep; }

}
}
#endif 
