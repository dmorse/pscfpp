#ifndef RP_MAX_ORDER_PARAMETER_TPP
#define RP_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   MaxOrderParameter<D,T>::MaxOrderParameter(
                                     Simulator<D,T>& simulator,
                                     System<D,T>& system)
    : AverageAnalyzer<D,T>(simulator, system),
      kSize_(-1)
   {  ParamComposite::setClassName("MaxOrderParameter"); }

   /*
   * Setup before main loop.
   */
   template <int D, class T>
   void MaxOrderParameter<D,T>::setup()
   {
      // Precondition: Require that the system has two monomer types
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      AverageAnalyzer<D,T>::setup();

      // Set mesh dimensions
      meshDimensions_ = system().domain().mesh().dimensions();
      FFTT::computeKMesh(meshDimensions_, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!wK_.isAllocated()){
         UTIL_CHECK(!psi_.isAllocated());
         wK_.allocate(meshDimensions_);
         psi_.allocate(kMeshDimensions_);
      }
      UTIL_CHECK(wK_.capacity() == kSize_);
      UTIL_CHECK(psi_.capacity() == kSize_);
   }

   /*
   * Compute array psi_ of squared Fourier amplitudes.
   */
   template <int D, class T>
   void MaxOrderParameter<D,T>::computePsi()
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
   template <int D, class T>
   void 
   MaxOrderParameter<D,T>::findMaximum(Array<typename T::Real> const & psi)
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
   template <int D, class T>
   void MaxOrderParameter<D,T>::outputValue(int step, double value)
   {
      std::ofstream& file = AverageAnalyzer<D,T>::outputFile_;
      UTIL_CHECK(file.is_open());
      int nSamplePerOutput = AverageAnalyzer<D,T>::nSamplePerOutput();
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
