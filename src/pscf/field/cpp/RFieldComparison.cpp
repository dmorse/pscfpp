/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.h"

// Explicit instantiation definitions
namespace Pscf {

   template class RFieldComparison<1,CPT>;
   template class RFieldComparison<2,CPT>;
   template class RFieldComparison<3,CPT>;

}
