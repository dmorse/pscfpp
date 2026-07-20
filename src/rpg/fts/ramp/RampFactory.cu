/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RampFactory.h"  

#include "LinearRamp.h"

#include <rp/fts/ramp/RampFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RampFactory<1, CudaTp<1> >;
      template class RampFactory<2, CudaTp<2> >;
      template class RampFactory<3, CudaTp<3> >;
   }
}
