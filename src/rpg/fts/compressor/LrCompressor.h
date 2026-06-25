#ifndef RPG_LR_COMPRESSOR_H
#define RPG_LR_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrCompressor.h>      // direct base template
#include <rpg/system/Types.h>                    // direct base argument
#include <rpg/fts/compressor/IntraCorrelation.h> // direct base member
#include <prdc/cuda/RField.h>                    // direct base member
#include <prdc/cuda/RFieldDft.h>                 // direct base member
#include <rp/fts/compressor/Compressor.h>        // indirect base class

namespace Pscf {
namespace Rpg {

   // Namespaces that can be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * Linear-response Anderson mixing compressor.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::LrCompressor, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::LrCompressor
   * \see \ref rp_LrCompressor_page "Manual Page"
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class LrCompressor
    : public Rp::LrCompressor<D, Rpg::Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LrCompressor(Rp::System<D, Rpg::Types<D> >& system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LrCompressor<1, Rpg::Types<1> >;
      extern template class LrCompressor<2, Rpg::Types<2> >;
      extern template class LrCompressor<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class LrCompressor<1>;
      extern template class LrCompressor<2>;
      extern template class LrCompressor<3>;
   }
}
#endif
