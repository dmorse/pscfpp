/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/FftwDRArray.h>
#include <pscf/cpu/send.h>
#include <pscf/backend/CPT.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1,CPT>;
      template class IntraCorrelation<2,CPT>;
      template class IntraCorrelation<3,CPT>;
   }
}
