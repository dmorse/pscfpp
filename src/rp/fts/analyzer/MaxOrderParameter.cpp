/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <pscf/backend/cpp/VecOpCx.h>

#include <rp/fts/analyzer/MaxOrderParameter.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   // Explicit instantiation definitions
   //PSCF_TMPL_DEFINE_CPP(MaxOrderParameter)
   template class MaxOrderParameter<1,CPT>;
   template class MaxOrderParameter<2,CPT>;
   template class MaxOrderParameter<3,CPT>;

} // namespace Rp
} // namespace Pscf
