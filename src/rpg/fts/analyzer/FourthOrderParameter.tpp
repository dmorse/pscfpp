#ifndef RPG_FOURTH_ORDER_PARAMETER_TPP
#define RPG_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/resources.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D>::FourthOrderParameter(Simulator<D>& simulator,
                                                 System<D>& system)
    : AverageAnalyzer<D>(simulator, system),
      kSize_(1),
      isInitialized_(false)
   {  setClassName("FourthOrderParameter"); }

   /*
   * Destructor.
   */
   template <int D>
   FourthOrderParameter<D>::~FourthOrderParameter()
   {}

   /*
   * FourthOrderParameter setup
   */
   template <int D>
   void FourthOrderParameter<D>::setup()
   {
      // Local copies of data
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().mesh().dimensions();

      // Precondition: Require that the system has two monomer types
      UTIL_CHECK(nMonomer == 2);

      AverageAnalyzer<D>::setup();

      // Compute DFT k-space mesh kMeshDimensions_ and kSize_
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wc0_.allocate(dimensions);
         wK_.allocate(dimensions);
         prefactor_.allocate(kMeshDimensions_);
         VecOp::eqS(prefactor_, 0);
      }

      isInitialized_ = true;

      computePrefactor();
   }

   template <int D>
   double FourthOrderParameter<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());

      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      RField<D> psi;
      psi.allocate(kMeshDimensions_);

      // Convert W_(r) to fourier mode W_(k)
      VecOp::eqV(wc0_, simulator().wc(0));
      system().fft().forwardTransform(wc0_, wK_);

      // psi = |wK_|^4
      VecOp::sqSqNormV(psi, wK_);

      // W_(k)^4 * weight factor
      VecOp::mulEqV(psi, prefactor_);

      // Get sum over all wavevectors
      FourthOrderParameter_ = Reduce::sum(psi);
      FourthOrderParameter_ = std::pow(FourthOrderParameter_, 0.25);

      return FourthOrderParameter_;
   }

   template <int D>
   void FourthOrderParameter<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double chi = system().interaction().chi(0,1);

         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(chi);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }

   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
   {
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      HostDArray<cudaReal> prefactor_h(kSize_);
      for (int i = 0; i < kSize_; ++i){
         prefactor_h[i] = 0;
      }
      IntVec<D> G;
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      DArray<IntVec<D>> GminList;
      GminList.allocate(kSize_);
      MeshIterator<D> itr(kMeshDimensions_);
      MeshIterator<D> searchItr(kMeshDimensions_);

      // Calculate GminList
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         GminList[itr.rank()] = Gmin;
      }

      // Compute weight factor for each G wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If the weight factor of the current wavevector has not been assigned
         if (prefactor_h[itr.rank()] == 0){
            Gmin = GminList[itr.rank()];

            // Compute inverse of wavevector
            nGmin.negate(Gmin);

            // Search for inverse of wavevector
            searchItr = itr;
            for (; !searchItr.atEnd(); ++searchItr){
               if (nGmin == GminList[searchItr.rank()]){
                  prefactor_h[itr.rank()] = 1.0/2.0;
                  prefactor_h[searchItr.rank()] = 1.0/2.0;
                  inverseFound = true;
               }
            }

            if (inverseFound == false){
               prefactor_h[itr.rank()]  = 1.0;
            }

         }

      }

      // Copy the weight factor from cpu(host) to gpu(device)
      prefactor_ = prefactor_h;
   }

}
}
#endif
