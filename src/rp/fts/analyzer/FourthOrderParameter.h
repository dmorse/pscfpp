#ifndef RP_FOURTH_ORDER_PARAMETER_H
#define RP_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h>   // base class template
#include <pscf/math/IntVec.h>                  // member

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * FourthOrderParameter is used to detect an order-disorder transition.
   *
   * This class evaluates and averages the sum of fourth power of the
   * magnitude of the Fourier mode amplitudes of a fluctuating exchange
   * field for a binary system.
   *
   * The order parameter is defined as
   * \f[
   *     \Psi_{\text{fourth}} \equiv
   *     \left[ \sum_{\bf G} W_{-}({\bf G})^4 \right] ^{\frac{1}{4}}
   * \f]
   * where \f$W_({\bf G})\f$ is a Fourier mode of fluctuating field.
   *
   * Specializations of this template are used as base classes for two
   * closely analogous class templates, also named FourthOrderParameter, 
   * that are defined in the Rpc and Rpg namespaces for use in the 
   * pscf_rpc and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_FourthOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class FourthOrderParameter : public T::AverageAnalyzer
   {

   public:

      /**
      * Setup before the main loop.
      */
      void setup() override;

   protected:

      /// Prefactor for each Fourier component.
      typename T::RField prefactor_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(typename T::Simulator& simulator, 
                           System<D,T>& system);

      /**
      * Destructor.
      */
      ~FourthOrderParameter() = default;

      /**
      * Compute and return the order parameter.
      */
      double compute() override;

      /**
      * Compute prefactor for each Fourier wavevector.
      *
      * For the real-valued function W_, each Fourier
      * coefficient G satisfies W_(G) = W_(-G). This function
      * uses Brillouin Zone (BZ) indices representation. After
      * applying fftw, if both the wavevector G and its
      * inverse -G exist in k-space, the prefactor is
      * assigned to be 1/2 for both G and -G. Otherwise,
      * it is assigned to be 1.
      */
      void computePrefactor(Array<double>& prefactor);

      // Alias for base class
      using AverageAnalyzerT = typename T::AverageAnalyzer;

      // Inherited protected member functions (selected).
      using AverageAnalyzerT::simulator;
      using AverageAnalyzerT::system;

   private:

      /// Fourier transform of W_ field.
      typename T::RFieldDft wK_;

      /// Fourth powers of Fourier magnitudes, with prefactors.
      typename T::RField psi_;

      /// Dimensions of Fourier space (k-grid) mesh for a real field.
      IntVec<D> kMeshDimensions_;

      /// Has setup been completed?
      bool isInitialized_;

      /**
      * Initialize prefactor_ member array.
      *
      * The GPU version of this function must compute values on
      * on the CPU host and then copy them to a device array. 
      */
      virtual void computePrefactor() = 0;

      using FFTT = typename T::FFT;

   };

}
}
#endif
