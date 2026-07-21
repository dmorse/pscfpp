#ifndef RP_INTRACORRELATION_H
#define RP_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>         // member
#include <util/containers/DArray.h>   // member

// Forward declarations
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class RField;
   }
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
   }
   namespace Correlation {
      template <typename WT> class Mixture;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Intramolecular correlation analyzer.
   *
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T>
   class IntraCorrelation
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      IntraCorrelation(System<D,T> const & system);

      /**
      * Destructor.
      */
      ~IntraCorrelation();

      /**
      * Compute total intramolecular correlation function (all blocks).
      *
      * \param correlations  omega values on a k-space mesh
      */
      virtual
      void computeOmegaTotal(RField<D,T>& correlations);

   protected:

      /**
      * Compute total intramolecular correlation function (all blocks).
      *
      * \param correlations  omega values on a k-space mesh
      */
      void computeOmegaTotalArray(Array<typename T::Real>& correlations);

      /**
      * Get and store r-grid and kgrid-mesh dimensions.
      *
      * Results are stored in private member variables for later use.
      */
      void getMeshDimensions();

      /**
      * Compute array of squared wavevectors.
      *
      * Results are stored in a private member variable.
      */
      void computeGsq();

      /*
      * Get the size of the k-space mesh. 
      */
      int kSize() const;

      /**
      * Get the parent system by const ref.
      */
      System<D,T> const & system() const;

   private:

      using RealT = typename T::Real;

      /// Pointer to parent system object.
      System<D,T> const * systemPtr_;

      /// Pointer to child Correlation::Mixture object.
      Correlation::Mixture<RealT>* correlationMixturePtr_;

      /// Array of square magnitudes for wavevectors on a k-grid.
      DArray<RealT> Gsq_;

      /// Local (host) array of omega values on a k-grid.
      typename T::RLocArray correlations_;

      /// Dimensions of r-space mesh.
      IntVec<D> meshDimensions_;

      /// Dimensions of k-space mesh for DFT of a real function.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in the k-space mesh.
      int kSize_;

      using FFTT = typename T::FFT;

   };

   // Get the number of wavevectors in the k-space mesh.
   template <int D, class T> inline
   int IntraCorrelation<D,T>::kSize() const
   {  return kSize_; }

   // Get the parent system by const reference.
   template <int D, class T> inline
   System<D,T> const & IntraCorrelation<D,T>::system() const
   {  return *systemPtr_; }

} // namespace Rp
} // namespace Pscf
#endif
