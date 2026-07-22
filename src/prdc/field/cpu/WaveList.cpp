/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.tpp"

// Explicit instantiation definitions
namespace Pscf { 
   namespace Prdc { 
      template class WaveList<1,CPT>;
      template class WaveList<2,CPT>;
      template class WaveList<3,CPT>;
   }
}
