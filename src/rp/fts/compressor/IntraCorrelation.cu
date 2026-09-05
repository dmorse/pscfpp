/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/send.h>
#include <pscf/backend/cuda/HostDArray.h>
#include <pscf/backend/cuda/CUT.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1,CUT>;
      template class IntraCorrelation<2,CUT>;
      template class IntraCorrelation<3,CUT>;
   }
}
