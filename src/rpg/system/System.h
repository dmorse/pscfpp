#ifndef RPG_SYSTEM_H
#define RPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <rp/system/System.h>       // base class template
#include <pscf/cuda/CudaTp.h>       // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::System<1, CudaTp<1> >;
      extern template class Rp::System<2, CudaTp<1> >;
      extern template class Rp::System<3, CudaTp<1> >;
   }
}
#endif
