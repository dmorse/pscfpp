/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Domain.tpp>
#include <pscf/backend/CPT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Domain<1,CPT>;
      template class Domain<2,CPT>;
      template class Domain<3,CPT>;
   }
}
