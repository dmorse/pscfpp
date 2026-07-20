/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/scft/sweep/SweepFactory.h>     // header
#include <rpg/scft/sweep/LinearSweep.h>
#include <rp/scft/sweep/SweepFactory.tpp>    // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepFactory<1, CudaTp<1> >;
      template class SweepFactory<2, CudaTp<2> >;
      template class SweepFactory<3, CudaTp<3> >;
   }
}
