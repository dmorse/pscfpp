#ifndef RPG_PERTURBATION_H
#define RPG_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <rpg/system/Types.h>                  // base class argument

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
      extern template class Perturbation<1, Rpg::Types<1> >;
      extern template class Perturbation<2, Rpg::Types<2> >;
      extern template class Perturbation<3, Rpg::Types<3> >;
   }
}
#endif
