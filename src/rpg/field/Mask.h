#ifndef RPG_MASK_H
#define RPG_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Mask.h>            // lass template
#include <pscf/cuda/Cuda.h>         // class template argument
#include <prdc/field/cuda/RField.h>   // base class member


// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Mask< 1, CudaTp<1> >;
      extern template class Mask< 2, CudaTp<2> >;
      extern template class Mask< 3, CudaTp<3> >;
   }
}
#endif
