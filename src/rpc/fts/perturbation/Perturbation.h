#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <rpc/system/Types.h>                  // base class argument

// Forward declarations
namespace Pscf {
   namespace Prdc {
      namespace Cpu {
         template <int D> class RField;
      }
   }
   namespace Rp {
      template <int D, class T> class System;
   }
   namespace Rpc {
      template <int D> class Simulator;
   }
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1, Rpc::Types<1> >;
      extern template class Perturbation<2, Rpc::Types<2> >;
      extern template class Perturbation<3, Rpc::Types<3> >;
   }
}
#endif
