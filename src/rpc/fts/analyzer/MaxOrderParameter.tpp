#ifndef RPC_MAX_ORDER_PARAMETER_TPP
#define RPC_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/cpu/VecOpCx.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   MaxOrderParameter<D>::MaxOrderParameter(Simulator<D>& simulator,
                                           System<D>& system)
    : AverageAnalyzer<D>(simulator, system),
      kSize_(1),
      isInitialized_(false)
   {  ParamComposite::setClassName("MaxOrderParameter"); }

   /*
   * Setup before main loop.
   */
   template <int D>
   void MaxOrderParameter<D>::setup()
   {

      // Precondition: Require that the system has two monomer types
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      AverageAnalyzer<D>::setup();

      // Set mesh dimensions
      meshDimensions_ = system().domain().mesh().dimensions();
      FFT<D>::computeKMesh(meshDimensions_, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wK_.allocate(meshDimensions_);
         psi_.allocate(kMeshDimensions_);
      }

      isInitialized_ = true;
   }

   /*
   * Search for and return maximum Fourier amplitude.
   */
   template <int D>
   double MaxOrderParameter<D>::compute()
   {
      computePsi();
      findMaximum(psi_);
      return maxPsi_;
   }

   /*
   * Compute array psi_ of squared Fourier amplitudes.
   */
   template <int D>
   void MaxOrderParameter<D>::computePsi()
   {
      UTIL_CHECK(system().w().hasData());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);
      VecOp::sqAbsV(psi_, wK_);
   }

   /*
   * Search for and return maximum Fourier amplitude.
   */
   template <int D>
   void MaxOrderParameter<D>::findMaximum(Array<double> const & psi)
   {
      // Identify index of maximum element of array psi
      maxPsi_ = psi[1];
      int maxIndex = 1;
      for (int i = 2; i < kSize_; ++i){
         if (psi[i] > maxPsi_){
            maxPsi_ = psi[i];
            maxIndex = i;
         }
      }

      // Find minimal indices of wavevector Gmax_ of maximum element
      Mesh<D> kMesh(kMeshDimensions_);
      Gmax_ = kMesh.position(maxIndex);
      UnitCell<D> const & unitCell = system().domain().unitCell();
      Gmax_ = shiftToMinimum(Gmax_, meshDimensions_, unitCell);
   }

   /*
   * Output instantaneous value during simulation.
   */
   template <int D>
   void MaxOrderParameter<D>::outputValue(int step, double value)
   {
      std::ofstream& file = AverageAnalyzer<D>::outputFile_;
      UTIL_CHECK(file.is_open());
      int nSamplePerOutput = AverageAnalyzer<D>::nSamplePerOutput();
      if (nSamplePerOutput == 1) {
         file << Int(step);
         file << "   ( ";
         for (int i = 0; i < D; i++){
            file << Int(Gmax_[i],3) << " ";
         }
         file << " )  ";
         file << Dbl(value);
         file << "\n";
      } else {
         file << Int(step);
         file << Dbl(value);
         file << "\n";
      }
   }

}
}
#endif
