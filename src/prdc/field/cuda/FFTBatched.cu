/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.tpp"

// Explicit instantiation definitions
namespace Pscf {
   namespace Prdc {
      template class FFTBatched<1>;
      template class FFTBatched<2>;
      template class FFTBatched<3>;
   } 
} 
