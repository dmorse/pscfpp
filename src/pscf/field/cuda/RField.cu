/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.tpp"

// Explicit instantiation definitions
namespace Pscf {
   template class RField<1,CUT>;
   template class RField<2,CUT>;
   template class RField<3,CUT>;
}
