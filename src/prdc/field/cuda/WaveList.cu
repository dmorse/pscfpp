/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.tpp"

// Explicit instantation definitions
namespace Pscf { 
   namespace Prdc { 
      template class WaveList<1,CUT>;
      template class WaveList<2,CUT>;
      template class WaveList<3,CUT>;
   }
}
