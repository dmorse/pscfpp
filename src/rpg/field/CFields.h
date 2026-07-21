#ifndef RPG_C_FIELDS_H
#define RPG_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>     // class template
#include <pscf/cuda/CudaTp.h>     // class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      extern template class Rp::CFields<1, CudaTp<1> >;
      extern template class Rp::CFields<2, CudaTp<2> >;
      extern template class Rp::CFields<3, CudaTp<3> >;
   }
}
#endif
