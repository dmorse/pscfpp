/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/cpu/RField.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/fts/compressor/AmCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AmCompressor<1,CPT>;
      template class AmCompressor<2,CPT>;
      template class AmCompressor<3,CPT>;
   }
}
