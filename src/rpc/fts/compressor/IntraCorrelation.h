#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>         // memmber variable type
#include <util/containers/DArray.h>   // memmber variable type



// Forward declarations
namespace Pscf {
   namespace Correlation {
      template <typename WT> class Mixture;
   }
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
   class IntraCorrelation
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(System<D> const & system);

      /**
      * Destructor.
      */
      ~IntraCorrelation();

      /**
      * Compute total intramolecular correlation function (all blocks).
      *
      * \param correlations  omega values on a k-space mesh
      */
      void computeOmegaTotal(RField<D>& correlations);

   protected:

      /**
      * Compute k-grid mesh dimensions and array of squared wavevectors.
      *
      * Results are stored in private member variables for later use.
      */
      void computeMeshProperties();

      /**
      * Get the parent system by const ref.
      */
      System<D> const & system() const;

   private:

      /// Pointer to parent system object.
      System<D> const * systemPtr_;

      /// Pointer to child Correlation::Mixture object.
      Correlation::Mixture<double>* correlationMixturePtr_;

      /// Array of square magnitudes for wavevectors on a k-grid.
      DArray<double> Gsq_;

      /// Dimensions of k-space mesh for DFT of a real function.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in the k-space mesh.
      int kSize_;

   };

   // Get the parent system by const reference.
   template <int D> inline
   System<D> const & IntraCorrelation<D>::system() const
   {  return *systemPtr_; }

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
#include <pscf/correlation/Mixture.h>
namespace Pscf {
   namespace Correlation {
      extern template class Mixture<double>;
   }
   namespace Rpc {
      extern template class IntraCorrelation<1>;
      extern template class IntraCorrelation<2>;
      extern template class IntraCorrelation<3>;
   }
}
#endif
