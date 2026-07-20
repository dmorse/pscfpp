/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearSweep.h"
#include <rpg/system/System.h>

#include <rp/scft/sweep/LinearSweep.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LinearSweep< 1, CudaTp<1> >;
      template class LinearSweep< 2, CudaTp<2> >;
      template class LinearSweep< 3, CudaTp<3> >;
   }
}
