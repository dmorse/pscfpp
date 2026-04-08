#ifndef RPG_LR_AM_COMPRESSOR_H
#define RPG_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrAmCompressor.h>    // direct base template
#include <rpg/system/Types.h>                    // direct base argument
#include <rpg/fts/compressor/AmCompressorBase.h> // indirect base class
#include <rpg/fts/compressor/IntraCorrelation.h> // direct base member
#include <prdc/cuda/RField.h>                    // direct base member
#include <prdc/cuda/RFieldDft.h>                 // direct base member

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Linear-response Anderson mixing compressor.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of the base class template Rp::LrAmCompressor, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::LrAmCompressor
   * \see \ref rp_LrAmCompressor_page "Manual Page"
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class LrAmCompressor
    : public Rp::LrAmCompressor<D, Rpg::Types<D>, DeviceArray<cudaReal> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LrAmCompressor(System<D>& system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class LrAmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      extern template 
      class LrAmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      extern template 
      class LrAmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
   namespace Rpg {
      extern template class LrAmCompressor<1>;
      extern template class LrAmCompressor<2>;
      extern template class LrAmCompressor<3>;
   }
}
#endif
