#ifndef RPC_RAMP_TPP
#define RPC_RAMP_TPP

#include "Ramp.h"
#include <rpc/simulate/Simulator.h>
#include <prdc/cpu/RField.h>

#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

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

   #if 0   
   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D>
   void Ramp<D>::readParameters(std::istream& in)
   {}
   #endif

}
}
#endif 
