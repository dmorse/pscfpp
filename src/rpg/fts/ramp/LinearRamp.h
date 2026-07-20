#ifndef RPG_LINEAR_RAMP_H
#define RPG_LINEAR_RAMP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/LinearRamp.h>      // direct base class template
#include <pscf/cuda/Cuda.h>            // base class template argument
#include <rpg/fts/ramp/RampParameter.h>  // base class member
#include <rpg/fts/ramp/Ramp.h>           // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LinearRamp<1, CudaTp<1> >;
      extern template class LinearRamp<2, CudaTp<2> >;
      extern template class LinearRamp<3, CudaTp<3> >;
   }
}
#endif
