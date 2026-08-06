/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/Reduce.h>
#include <pscf/backends/CUT.h>

#include <rp/scft/ScftThermo.tpp>      // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1,CUT>;
      template class ScftThermo<2,CUT>;
      template class ScftThermo<3,CUT>;
   }
}
