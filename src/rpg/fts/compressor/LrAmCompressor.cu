/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <prdc/field/cuda/FFT.h>
#include <prdc/field/cuda/RField.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/VecOpMisc.h>
#include <pscf/cuda/Reduce.h>

#include <rp/fts/compressor/LrAmCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrAmCompressor<1,CUT>;
      template class LrAmCompressor<2,CUT>;
      template class LrAmCompressor<3,CUT>;
   }
}
