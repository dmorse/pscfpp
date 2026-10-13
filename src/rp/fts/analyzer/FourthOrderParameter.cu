/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/FourthOrderParameter.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FourthOrderParameter<1,CUT>;
      template class FourthOrderParameter<2,CUT>;
      template class FourthOrderParameter<3,CUT>;
   }
}
