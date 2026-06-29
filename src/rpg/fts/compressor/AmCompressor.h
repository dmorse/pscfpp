#ifndef RPG_AM_COMPRESSOR_H
#define RPG_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/AmCompressor.h>      // direct base template
#include <rpg/system/Types.h>                    // direct base argument
#include <rpg/fts/compressor/AmCompressorBase.h> // indirect base 
#include <prdc/field/cuda/RField.h>                    // direct base member

namespace Pscf {
namespace Rpg {

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Anderson mixing compressor.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::AmCompressor, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::AmCompressor
   * \see \ref rp_AmCompressor_page "Manual Page"
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor
    : public Rp::AmCompressor<D, Rpg::Types<D>, DeviceArray<cudaReal> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      AmCompressor(Rp::System<D, Rpg::Types<D> >& system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class AmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      extern template 
      class AmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      extern template 
      class AmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
   namespace Rpg {
      extern template class AmCompressor<1>;
      extern template class AmCompressor<2>;
      extern template class AmCompressor<3>;
   }
}
#endif
