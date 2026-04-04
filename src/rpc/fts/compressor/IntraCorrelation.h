#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/IntraCorrelation.h>
#include <rpc/system/Types.h>

//#include <pscf/math/IntVec.h>         // memmber variable type
//#include <util/containers/DArray.h>   // memmber variable type

// Forward declarations
namespace Pscf {
   namespace Prdc {
      namespace Cpu {
         template <int D> class RField;
      }
   }
   namespace Rpc {
      template <int D> class System;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Intramolecular correlation analyzer.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class IntraCorrelation : public Rp::IntraCorrelation<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(System<D> const & system);

      #if 0
      /**
      * Compute total intramolecular correlation function (all blocks).
      *
      * \param correlations  omega values on a k-space mesh
      */
      void computeOmegaTotal(RField<D>& correlations);
      #endif

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class IntraCorrelation<1, Rpc::Types<1> >;
      extern template class IntraCorrelation<2, Rpc::Types<2> >;
      extern template class IntraCorrelation<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class IntraCorrelation<1>;
      extern template class IntraCorrelation<2>;
      extern template class IntraCorrelation<3>;
   }
}
#endif
