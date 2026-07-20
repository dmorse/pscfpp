#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <pscf/cpu/CppTp.h>                  // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1, CppTp<1> >;
      extern template class Perturbation<2, CppTp<2> >;
      extern template class Perturbation<3, CppTp<3> >;
   }
}
#endif
