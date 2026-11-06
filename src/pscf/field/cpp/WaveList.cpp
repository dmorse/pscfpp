/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.tpp"

// Explicit instantiation definitions
namespace Pscf { 
   template class WaveList<1,CPT>;
   template class WaveList<2,CPT>;
   template class WaveList<3,CPT>;
}
