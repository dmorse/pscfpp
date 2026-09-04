/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.tpp"           // class template implementation
#include <pscf/backend/CPT.h>   // backend type class

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class CFields<1,CPT>;
      template class CFields<2,CPT>;
      template class CFields<3,CPT>;
   }
}
