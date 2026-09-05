/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CudaVecRandom.h>
#include <pscf/backend/cuda/CUT.h>

#include <rp/fts/montecarlo/McSimulator.tpp>  // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McSimulator<1,CUT>;
      template class McSimulator<2,CUT>;
      template class McSimulator<3,CUT>;
   }
}
