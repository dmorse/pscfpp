#ifndef RPC_AM_COMPRESSOR_H
#define RPC_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/AmCompressor.h>      // direct base template
#include <rpc/system/Types.h>                    // direct base argument
#include <prdc/field/cpu/RField.h>                     // direct base member
#include <rpc/fts/compressor/AmCompressorBase.h> // indirect base

namespace Pscf {
namespace Rpc {

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Anderson mixing compressor.
   *
   * \see \ref rp_AmCompressor_page "Manual Page"
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor
    : public Rp::AmCompressor<D, Rpc::Types<D>, DArray<double> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      AmCompressor(Rp::System<D, Rpc::Types<D> >& system);

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AmCompressor<1, Rpc::Types<1>, DArray<double> >;
      extern template class AmCompressor<2, Rpc::Types<2>, DArray<double> >;
      extern template class AmCompressor<3, Rpc::Types<3>, DArray<double> >;
   }
   namespace Rpc {
      extern template class AmCompressor<1>;
      extern template class AmCompressor<2>;
      extern template class AmCompressor<3>;
   }
}
#endif
