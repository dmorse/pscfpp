#ifndef RPG_PERTURBATION_H
#define RPG_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <pscf/cuda/Cuda.h>                  // base class argument

// Forward declarations
namespace Pscf {
   namespace Prdc {
      namespace Cuda {
         template <int D> class RField;
      }
   }
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1, CudaTp<1> >;
      extern template class Perturbation<2, CudaTp<2> >;
      extern template class Perturbation<3, CudaTp<3> >;
   }
}
#endif
