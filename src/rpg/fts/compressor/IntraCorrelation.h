#ifndef RPG_INTRACORRELATION_H
#define RPG_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/IntraCorrelation.h>
#include <pscf/backends/CUT.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaTypes.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class IntraCorrelation<1,CUT>;
      extern template class IntraCorrelation<2,CUT>;
      extern template class IntraCorrelation<3,CUT>;
   }
}
#endif
