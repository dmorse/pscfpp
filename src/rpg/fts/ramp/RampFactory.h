#ifndef RPC_RAMP_FACTORY_H
#define RPC_RAMP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/RampFactory.h>
#include <pscf/cuda/Cuda.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RampFactory<1, CudaTp<1> >;
      extern template class RampFactory<2, CudaTp<2> >;
      extern template class RampFactory<3, CudaTp<3> >;
   }
}
#endif
