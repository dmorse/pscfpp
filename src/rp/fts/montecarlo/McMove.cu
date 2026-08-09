/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>
#include <rp/fts/montecarlo/McMove.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMove<1,CUT>;
      template class McMove<2,CUT>;
      template class McMove<3,CUT>;
   }
}
