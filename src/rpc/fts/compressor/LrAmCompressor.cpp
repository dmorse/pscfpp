/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>

#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <rp/fts/compressor/LrAmCompressor.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrAmCompressor<1, CppTp<1> >;
      template class LrAmCompressor<2, CppTp<2> >;
      template class LrAmCompressor<3, CppTp<3> >;
   }
}
