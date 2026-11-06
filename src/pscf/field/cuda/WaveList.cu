/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.tpp"

// Explicit instantation definitions
namespace Pscf { 
   template class WaveList<1,CUT>;
   template class WaveList<2,CUT>;
   template class WaveList<3,CUT>;
}
