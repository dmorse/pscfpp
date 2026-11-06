/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/MaxOrderParameter.tpp>
#include <pscf/backend/cpp/CPT.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Explicit instantiation definitions
   //PSCF_TMPL_DEFINE_CPP(MaxOrderParameter)
   template class MaxOrderParameter<1,CPT>;
   template class MaxOrderParameter<2,CPT>;
   template class MaxOrderParameter<3,CPT>;

} // namespace Rp
} // namespace Pscf
