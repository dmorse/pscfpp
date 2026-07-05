#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <rpc/system/Types.h>                  // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1, Rpc::Types<1> >;
      extern template class Perturbation<2, Rpc::Types<2> >;
      extern template class Perturbation<3, Rpc::Types<3> >;
   }
}
#endif
