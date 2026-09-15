#ifndef RP_FOURTH_ORDER_PARAMETER_H
#define RP_FOURTH_ORDER_PARAMETER_H

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

#include <pscf/backend/TmplDeclare.h>

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
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_FourthOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class FourthOrderParameter : public AverageAnalyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(Simulator<D,T>& simulator, 
                           System<D,T>& system);

      /**
      * Destructor.
      */
      ~FourthOrderParameter() = default;

      /**
      * Setup before the main loop.
      */
      void setup() override;

   protected:

      /**
      * Compute and return the order parameter.
      */
      double compute() override;

      // Inherited protected member functions (selected).
      using AverageAnalyzer<D,T>::simulator;
      using AverageAnalyzer<D,T>::system;

   private:

      /// Fourier transform of W_ field.
      RFieldDft<D,T> wK_;

      /// Fourth powers of Fourier magnitudes, with prefactors.
      RField<D,T> psi_;

      /// Prefactor for each Fourier component.
      RField<D,T> prefactor_;

      /// Dimensions of Fourier space (k-grid) mesh for a real field.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      /// Has setup been completed?
      bool isInitialized_;

      /**
      * Initialize member variable prefactor_ (device array).
      */
      void computePrefactor();

      /**
      * Compute prefactor for each Fourier wavevector (on host).
      *
      * For the real-valued function W_, each Fourier
      * coefficient G satisfies W_(G) = W_(-G). This function
      * uses Brillouin Zone (BZ) index representation. 
      * If both the wavevector G and its inverse -G are in the
      * half-spaced used for a real field, the prefactor is
      * assigned to be 1/2 for both G and -G. Otherwise, the
      * prefactor assigned to be 1.
      */
      void computePrefactor(Array<double>& prefactor);

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(FourthOrderParameter)

}
}
#endif
