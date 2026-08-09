/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>
#include <rp/fts/ramp/RampParameter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RampParameter<1,CUT>;
      template class RampParameter<2,CUT>;
      template class RampParameter<3,CUT>;
   }
}
