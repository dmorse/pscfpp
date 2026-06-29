#ifndef RPC_LR_AM_COMPRESSOR_H
#define RPC_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrAmCompressor.h>    // direct base template
#include <rpc/system/Types.h>                    // direct base argument
#include <rpc/fts/compressor/AmCompressorBase.h> // indirect base
#include <rpc/fts/compressor/IntraCorrelation.h> // direct base member
#include <prdc/field/cpu/RField.h>                     // direct base member
#include <prdc/field/cpu/RFieldDft.h>                  // direct base member

namespace Pscf {
namespace Rpc {

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Linear-response Anderson mixing compressor.
   *
   * \see \ref rp_LrAmCompressor_page "Manual Page"
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class LrAmCompressor
    : public Rp::LrAmCompressor<D, Rpc::Types<D>, DArray<double> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LrAmCompressor(Rp::System<D, Rpc::Types<D> >& system);

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class LrAmCompressor<1, Rpc::Types<1>, DArray<double> >;
      extern template 
      class LrAmCompressor<2, Rpc::Types<2>, DArray<double> >;
      extern template 
      class LrAmCompressor<3, Rpc::Types<3>, DArray<double> >;
   }
   namespace Rpc {
      extern template class LrAmCompressor<1>;
      extern template class LrAmCompressor<2>;
      extern template class LrAmCompressor<3>;
   }
}
#endif
