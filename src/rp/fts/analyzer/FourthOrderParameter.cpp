/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <pscf/backend/cpp/VecOpCx.h>
#include <pscf/backend/cpp/ReduceCx.h>

#include <rp/fts/analyzer/FourthOrderParameter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FourthOrderParameter<1,CPT>;
      template class FourthOrderParameter<2,CPT>;
      template class FourthOrderParameter<3,CPT>;
   }
}
