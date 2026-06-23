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
#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/HostDArray.h>

// Forward declarations
namespace Pscf {
   namespace Prdc {
      namespace Cuda {
         template <int D> class RField;
      }
   }
   namespace Rpg {
      template <int D> class System;
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
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from corresponding specializations of the base class template 
   * Rp::IntraCorrelation, and inherit most of their source code from
   * this base class.  
   *
   * \see Rp::IntraCorrelation
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class IntraCorrelation : public Rp::IntraCorrelation< D, Types<D> >
   {

   public:

      // Alias for base class
      using RpIntraCorrelation = Rp::IntraCorrelation<D,Types<D> >;

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(Rp::System<D, Rpg::Types<D> > const & system);

      /**
      * Compute total intramolecular correlation function (all blocks).
      *
      * \param correlations  omega values on a k-space mesh
      */
      void computeOmegaTotal(RField<D>& correlations);

      // Using declaration to avoid hiding name from base class.
      using RpIntraCorrelation::computeOmegaTotal;

   private:

      /// Host array of square magnitudes for wavevectors on a k-grid.
      HostDArray<cudaReal> correlations_;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class IntraCorrelation<1, Rpg::Types<1> >;
      extern template class IntraCorrelation<2, Rpg::Types<2> >;
      extern template class IntraCorrelation<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class IntraCorrelation<1>;
      extern template class IntraCorrelation<2>;
      extern template class IntraCorrelation<3>;
   }
}
#endif
