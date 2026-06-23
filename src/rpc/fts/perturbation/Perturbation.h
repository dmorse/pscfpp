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

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Additive perturbation of standard FTS Hamiltonian (base class).
   *
   * \see \ref psfts_perturb_page "Manual Page"
   * \ingroup Rpc_Fts_Perturbation_Module
   */
   template <int D>
   class Perturbation : public Rp::Perturbation<D, Rpc::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      Perturbation(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~Perturbation() = default;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1, Rpc::Types<1> >;
      extern template class Perturbation<2, Rpc::Types<2> >;
      extern template class Perturbation<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Perturbation<1>;
      extern template class Perturbation<2>;
      extern template class Perturbation<3>;
   }
}
#endif
