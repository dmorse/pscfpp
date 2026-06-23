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

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Additive perturbation of standard FTS Hamiltonian (base class).
   *
   * \see \ref psfts_perturb_page "Manual Page"
   * \ingroup Rpg_Fts_Perturbation_Module
   */
   template <int D>
   class Perturbation : public Rp::Perturbation<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      Perturbation(Rp::Simulator<D, Rpg::Types<D> >& simulator);

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
      extern template class Perturbation<1, Rpg::Types<1> >;
      extern template class Perturbation<2, Rpg::Types<2> >;
      extern template class Perturbation<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Perturbation<1>;
      extern template class Perturbation<2>;
      extern template class Perturbation<3>;
   }
}
#endif
