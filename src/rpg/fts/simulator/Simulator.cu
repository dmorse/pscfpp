/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/backends/CUT.h>  

#include <rp/fts/simulator/Simulator.tpp> 

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Simulator<1,CUT>;
      template class Simulator<2,CUT>;
      template class Simulator<3,CUT>;
   }
}
