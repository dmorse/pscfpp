/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "CFieldComparison.tpp"

// Explicit instantiation definitions
namespace Pscf {
   template class CFieldComparison<1,CPT>;
   template class CFieldComparison<2,CPT>;
   template class CFieldComparison<3,CPT>;
}
