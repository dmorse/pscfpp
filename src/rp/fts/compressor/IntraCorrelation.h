#ifndef RP_INTRACORRELATION_H
#define RP_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>         // memmber
#include <util/containers/DArray.h>   // memmber

// Forward declaration
namespace Pscf {
   namespace Correlation {
      template <typename WT> class Mixture;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Intramolecular correlation analyzer.
   *
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T>
   class IntraCorrelation
   {

   public:

      using RealT = typename T::Real;

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(typename T::System const & system);

      /**
      * Destructor.
      */
      ~IntraCorrelation();

      /**
      * Compute total intramolecular correlation function (all blocks).
      *
      * \param correlations  omega values on a k-space mesh
      */
      void computeOmegaTotal(Array<RealT>& correlations);

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
      typename T::System const & system() const;

   private:

      /// Pointer to parent system object.
      typename T::System const * systemPtr_;

      /// Pointer to child Correlation::Mixture object.
      Correlation::Mixture<RealT>* correlationMixturePtr_;

      /// Array of square magnitudes for wavevectors on a k-grid.
      DArray<RealT> Gsq_;

      /// Dimensions of k-space mesh for DFT of a real function.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in the k-space mesh.
      int kSize_;

      using FFTT = typename T::FFT;

   };

   // Get the parent system by const reference.
   template <int D, class T> inline
   typename T::System const & IntraCorrelation<D,T>::system() const
   {  return *systemPtr_; }

} // namespace Rp
} // namespace Pscf
#endif
