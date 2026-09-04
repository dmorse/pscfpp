#ifndef RP_MAX_ORDER_PARAMETER_BASE_H
#define RP_MAX_ORDER_PARAMETER_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"              // base class template
#include <prdc/field/RField.h>            // member
#include <prdc/field/RFieldDft.h>         // member
#include <pscf/math/IntVec.h>             // member
#include <util/containers/DArray.h>       // member

#include <pscf/backend/TmplDeclare.h>
#include <iostream>

// Forward declaration
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
   * Evaluate max of square magnitude of Fourier modes for exchange field.
   *
   * This class evaluates an average for the maximum of the squared Fourier 
   * mode amplitude for the fluctuating exchange field W_{-}(r) for an AB
   * system (nMonomer == 2).
   *
   * The quantity of interest is defined as
   * \f[
   *     \Psi  = \max [ |W_{-}({\bf k})|^{2} ]
   * \f]
   * where \f$ W_{-}({\bf k})\f$ is fluctuating field component with
   * wavevector \f$ {\bf k} \f$.
   *
   * Template parameters:
   *
   *   - D : dimension of space (D=1, 2, or 3)
   *   - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_MaxOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class MaxOrderParameterBase : public AverageAnalyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      MaxOrderParameterBase(Simulator<D,T>& simulator, 
                            System<D,T>& system);

      /**
      * Destructor.
      */
      ~MaxOrderParameterBase() = default;

      /**
      * Setup before simulation loop.
      */
      void setup() override;

   protected:

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      void outputValue(int step, double value) override;

      /**
      * Compute the psi_ array of squared Fourier coefficients.
      */
      void computePsi();

      /**
      * Find the wavevector of maximum Fourier magnitude.
      *
      * Results for the maximum square magnitude and the indices
      * of the wavevector for which this occured are stored in
      * maxPsi_ and Gmax_, respectively.
      *
      * \param psi  array of squared Fourier coefficients
      */
      void findMaximum(Array<typename T::Real> const & psi);

      /// Square magnitude |W_|^2 in Fourier space.
      RField<D,T> psi_;

      /// Maximum square magnitude (value of maximum element of psi_).
      double maxPsi_;

      /// Indices of wavevector with maximum magnitude.
      IntVec<D> Gmax_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      // Inherited protected member variables (selected).
      using AverageAnalyzer<D,T>::simulator;
      using AverageAnalyzer<D,T>::system;

   private:

      /// Fourier transform of W_ field
      RFieldDft<D,T> wK_;

      /// Dimensions of real space (r-grid)  mesh
      IntVec<D> meshDimensions_;

      /// Dimensions of Fourier space (k-grid) mesh for a real field
      IntVec<D> kMeshDimensions_;

      using FFTT = FFT<D,T>;
   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(MaxOrderParameterBase)

   // Primary template for subclasses
   template <int D, class T> class MaxOrderParameter;
}
}
#endif
