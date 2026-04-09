#ifndef RP_MAX_ORDER_PARAMETER_H
#define RP_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                      // base class
#include <pscf/math/IntVec.h>                     // member
#include <util/containers/DArray.h>               // member
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

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
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named MaxOrderParameter,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension of space (D=1, 2, or 3)
   *   - T : Types class (Rpc::Types<D> or Rpg::Types<D>)
   *
   * \see \ref rp_MaxOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class MaxOrderParameter : public T::AverageAnalyzer
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(typename T::Simulator& simulator, 
                        typename T::System& system);

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

      using AverageAnalyzerT = typename T::AverageAnalyzer;
      using AverageAnalyzerT::simulator;
      using AverageAnalyzerT::system;

   protected:

      /// Square magnitude |W_|^2 in Fourier space.
      typename T::RField psi_;

      /// Maximum square magnitude (value of maximum element of psi_).
      double maxPsi_;

      /// Indices of wavevector with maximum magnitude.
      IntVec<D> Gmax_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      /// Real type used in field containers.
      using RealT = typename T::Real;

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
      void findMaximum(Array<RealT> const & psi);

   private:

      /// Fourier transform of W_ field
      typename T::RFieldDft wK_;

      /// Dimensions of real space (r-grid)  mesh
      IntVec<D> meshDimensions_;

      /// Dimensions of Fourier space (k-grid) mesh for a real field
      IntVec<D> kMeshDimensions_;

      using FFTT = typename T::FFT;
   };

}
}
#endif
