/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/MaxOrderParameter.tpp>
#include <pscf/backend/cuda/CUT.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Explicit instantiation definitions
   //PSCF_TMPL_DEFINE_CUDA(MaxOrderParameter)
   template class MaxOrderParameter<1,CUT>;
   template class MaxOrderParameter<2,CUT>;
   template class MaxOrderParameter<3,CUT>;

}
}

