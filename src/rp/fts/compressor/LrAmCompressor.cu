/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/IntraCorrelation.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <prdc/field/FFT.h>
#include <prdc/field/RField.h>

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/VecOpMisc.h>
#include <pscf/backend/cuda/Reduce.h>
#include <pscf/backend/cuda/CUT.h>

#include <rp/fts/compressor/LrAmCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrAmCompressor<1,CUT>;
      template class LrAmCompressor<2,CUT>;
      template class LrAmCompressor<3,CUT>;
   }
}
