/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/System.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/FFT.h>
#include <prdc/field/RField.h>

#include <pscf/math/IntVec.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/backends/CPT.h>

#include <util/global.h>

#include <rp/fts/compressor/IntraCorrelation.h>
#include <rp/fts/compressor/LrAmCompressor.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrAmCompressor<1,CPT>;
      template class LrAmCompressor<2,CPT>;
      template class LrAmCompressor<3,CPT>;
   }
}
