#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <pscf/backends/CPT.h>                  // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1,CPT>;
      extern template class Perturbation<2,CPT>;
      extern template class Perturbation<3,CPT>;
   }
}
#endif
