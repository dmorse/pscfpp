/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "MeshIteratorFortran.tpp"

// Explicit instantiation definitions
namespace Pscf {
   template class MeshIteratorFortran<1>;
   template class MeshIteratorFortran<2>;
   template class MeshIteratorFortran<3>;
}
