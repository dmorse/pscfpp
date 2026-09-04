/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/Reduce.h>
#include <pscf/backend/CPT.h>

#include <rp/field/Mask.tpp>        // class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::Mask<1,CPT>;
      template class Rp::Mask<2,CPT>;
      template class Rp::Mask<3,CPT>;
   }
}
