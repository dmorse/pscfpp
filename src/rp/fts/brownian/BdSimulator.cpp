/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/backends/CPT.h>

#include <rp/fts/brownian/BdSimulator.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdSimulator<1,CPT>;
      template class BdSimulator<2,CPT>;
      template class BdSimulator<3,CPT>;
   }
}
