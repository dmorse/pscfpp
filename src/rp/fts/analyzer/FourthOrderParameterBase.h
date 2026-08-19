#ifndef RP_FOURTH_ORDER_PARAMETER_BASE_H
#define RP_FOURTH_ORDER_PARAMETER_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h>   // base class template
#include <prdc/field/RField.h>                 // member
#include <prdc/field/RFieldDft.h>              // member
#include <pscf/math/IntVec.h>                  // member

#include <pscf/backends/TmplDeclare.h>

// Forward declarations
namespace Util {
   template <typename T> class Array;
}
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class FFT;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * FourthOrderParameterBase is used to detect an order-disorder transition.
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
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_FourthOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class FourthOrderParameterBase : public AverageAnalyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameterBase(Simulator<D,T>& simulator, 
                           System<D,T>& system);

      /**
      * Destructor.
      */
      ~FourthOrderParameterBase() = default;

      /**
      * Setup before the main loop.
      */
      void setup() override;

   protected:

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

      /// Prefactor for each Fourier component.
      RField<D,T> prefactor_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      // Alias for base class
      using AverageAnalyzerT = AverageAnalyzer<D,T>;

      // Inherited protected member functions (selected).
      using AverageAnalyzer<D,T>::simulator;
      using AverageAnalyzer<D,T>::system;

   private:

      /// Fourier transform of W_ field.
      RFieldDft<D,T> wK_;

      /// Fourth powers of Fourier magnitudes, with prefactors.
      RField<D,T> psi_;

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

      using FFTT = FFT<D,T>;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(FourthOrderParameterBase)

   // Primary template declaration for subclasses
   template <int D, class t> class FourthOrderParameter;

}
}
#endif
