#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/IntraCorrelation.h>   // base class template
#include <rpc/system/Types.h>                     // base class argument
#include <prdc/field/cpu/FftwDRArray.h>           // base class member

#if 0
namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Intramolecular correlation analyzer.
   *
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D>
   class IntraCorrelation : public Rp::IntraCorrelation<D, Rpc::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(Rp::System<D, Rpc::Types<D> > const & system);

   };

} // namespace Rpc
} // namespace Pscf
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class IntraCorrelation<1, Rpc::Types<1> >;
      extern template class IntraCorrelation<2, Rpc::Types<2> >;
      extern template class IntraCorrelation<3, Rpc::Types<3> >;
   }
}
#endif
