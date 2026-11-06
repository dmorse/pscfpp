/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.tpp"

// Explicit instantiation definitions
namespace Pscf {
   template class FFTBatched<1>;
   template class FFTBatched<2>;
   template class FFTBatched<3>;
} 
