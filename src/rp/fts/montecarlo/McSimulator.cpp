/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CpuVecRandom.h>
#include <pscf/backend/CPT.h>

#include <rp/fts/montecarlo/McSimulator.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McSimulator<1,CPT>;
      template class McSimulator<2,CPT>;
      template class McSimulator<3,CPT>;
   }
}
