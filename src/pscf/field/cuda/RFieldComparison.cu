/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.tpp"

// Explicit instantiation definitions
namespace Pscf {
   template class RFieldComparison<1,CUT>;
   template class RFieldComparison<2,CUT>;
   template class RFieldComparison<3,CUT>;
}
