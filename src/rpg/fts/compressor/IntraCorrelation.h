#ifndef RPG_INTRACORRELATION_H
#define RPG_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/IntraCorrelation.h>
#include <rpg/system/Types.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaTypes.h>

#if 0
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
   * Intramolecular correlation analyzer.
   *
   * \see Rp::IntraCorrelation
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class IntraCorrelation : public Rp::IntraCorrelation< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(Rp::System<D, Rpg::Types<D> > const & system);

   };

} // namespace Rpg
} // namespace Pscf
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class IntraCorrelation<1, Rpg::Types<1> >;
      extern template class IntraCorrelation<2, Rpg::Types<2> >;
      extern template class IntraCorrelation<3, Rpg::Types<3> >;
   }
}
#endif
