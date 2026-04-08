#ifndef RPC_MAX_ORDER_PARAMETER_H
#define RPC_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                      // base class
#include <prdc/cpu/RField.h>                      // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <util/containers/DArray.h>               // member
#include <iostream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * MaxOrderParameter is used to detect an order-disorder transition.
   *
   * This class evalutaes maximum of the squared Fourier mode amplitude
   * for the fluctuating field.
   *
   * The order parameter is defined as
   * \f[
   *     \psi(k)  = \max [ |W_{-}({\bf k})|^{2} ]
   * \f]
   * where \f$ W_{-}({\bf k})\f$ is fluctuating field component with
   * wavevector \f$ {\bf k} \f$.
   *
   * \see \ref rp_MaxOrderParameter_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class MaxOrderParameter : public AverageAnalyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**
      * Setup before simulation loop.
      */
      void setup() override;

   protected:

      /**
      * Compute and return the max order parameter.
      */
      double compute() override;

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      void outputValue(int step, double value) override;

      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;

   protected:

      /// Square magnitude |W_|^2 in Fourier space.
      RField<D> psi_;

      /// Maximum square magnitude (value of maximum element of psi_).
      double maxPsi_;

      /// Indices of wavevector with maximum magnitude.
      IntVec<D> Gmax_;

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
      void findMaximum(Array<double> const & psi);

   private:

      /// Fourier transform of W_ field
      RFieldDft<D> wK_;

      /// Dimensions of r-grid mesh
      IntVec<D> meshDimensions_;

      /// Dimensions of Fourier space (k-grid) mesh for a real field
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in Fourier space (k-grid) mesh.
      int  kSize_;

      /// Has setup been completed?
      bool isInitialized_;

   };

   // Explicit instantiation declarations
   extern template class MaxOrderParameter<1>;
   extern template class MaxOrderParameter<2>;
   extern template class MaxOrderParameter<3>;

}
}
#endif
