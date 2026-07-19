#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <pscf/cpu/Cpp.h>                  // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1, Cpp<1> >;
      extern template class Perturbation<2, Cpp<2> >;
      extern template class Perturbation<3, Cpp<3> >;
   }
}
#endif
