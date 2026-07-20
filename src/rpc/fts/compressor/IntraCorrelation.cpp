/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/send.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1, CppTp<1> >;
      template class IntraCorrelation<2, CppTp<2> >;
      template class IntraCorrelation<3, CppTp<3> >;
   }
}
