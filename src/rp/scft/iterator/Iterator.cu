/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>
#include <rp/scft/iterator/Iterator.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Iterator<1,CUT>;
      template class Iterator<2,CUT>;
      template class Iterator<3,CUT>;
   }
}
