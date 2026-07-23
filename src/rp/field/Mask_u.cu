/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <pscf/cuda/Reduce.h>
#include <rp/field/Mask.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Mask<1,CUT>;
      template class Mask<2,CUT>;
      template class Mask<3,CUT>;
   }
}
