/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.tpp"

// Explicit instantiation definitions
namespace Pscf {
   template class RFieldDft<1,CUT>;
   template class RFieldDft<2,CUT>;
   template class RFieldDft<3,CUT>;
}
