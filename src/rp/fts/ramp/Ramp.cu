/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>
#include <rp/fts/ramp/Ramp.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Ramp<1,CUT>;
      template class Ramp<2,CUT>;
      template class Ramp<3,CUT>;
   }
}
