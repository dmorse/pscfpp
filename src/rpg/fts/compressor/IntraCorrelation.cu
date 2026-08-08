/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>

#include <prdc/field/cuda/FFT.h>
#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/send.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1,CUT>;
      template class IntraCorrelation<2,CUT>;
      template class IntraCorrelation<3,CUT>;
   }
}
