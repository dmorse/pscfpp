/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrCompressor.h"
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <prdc/field/cpu/FFT.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <rp/fts/compressor/LrCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrCompressor<1,CPT>;
      template class LrCompressor<2,CPT>;
      template class LrCompressor<3,CPT>;
   }
}
