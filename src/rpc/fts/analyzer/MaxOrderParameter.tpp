#ifndef RPC_MAX_ORDER_PARAMETER_TPP
#define RPC_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/inter/Interaction.h>
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
   {  setClassName("MaxOrderParameter"); }

   /*
   * Destructor.
   */
   template <int D>
   MaxOrderParameter<D>::~MaxOrderParameter()
   {}

   /*
   * MaxOrderParameter setup
   */
   template <int D>
   void MaxOrderParameter<D>::setup()
   {

      // Precondition: Require that the system has two monomer types
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      AverageAnalyzer<D>::setup();

      // Set mesh dimensions
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wK_.allocate(dimensions);
      }

      isInitialized_ = true;

      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }

   /*
   * Search for and return maximum Fourier amplitude.
   */
   template <int D>
   double MaxOrderParameter<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());

      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      MeshIterator<D> itr;
      itr.setDimensions(kMeshDimensions_);
      std::vector<double> psi(kSize_);
      DArray<IntVec<D>> GminList;
      GminList.allocate(kSize_);
      IntVec<D> G;
      IntVec<D> Gmin;
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Conver W_(r) to fourier mode W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);

      for (itr.begin(); !itr.atEnd(); ++itr) {
         G = itr.position();
         Gmin = shiftToMinimum(G, dimensions, unitCell);
         GminList[itr.rank()] = Gmin;

         std::complex<double> wK(wK_[itr.rank()][0], wK_[itr.rank()][1]);
         psi[itr.rank()] = std::norm(wK);
      }

      maxOrderParameter_ = psi[1];
      int maxIndex = 1;
      for (int i = 2; i < kSize_; ++i){
         if (psi[i] > maxOrderParameter_){
            maxOrderParameter_ = psi[i];
            maxIndex = i;
         }
      }

      GminStar_ = GminList[maxIndex];

      return maxOrderParameter_;
   }

   /*
   * Output instantaneous value during simulation.
   */
   template <int D>
   void MaxOrderParameter<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double chi= system().interaction().chi(0,1);

         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(chi);
         outputFile_ << "   ( ";
         for (int i = 0; i < D; i++){
            outputFile_ << Int(GminStar_[i],3) << " ";
         }
         outputFile_ << " )  ";
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
      } else {
         //AverageAnalyzer<D>::outputValue(step, value);
         outputFile_ << Int(step);
         outputFile_ << "   ( ";
         for (int i = 0; i < D; i++){
            outputFile_ << Int(GminStar_[i],3) << " ";
         }
         outputFile_ << " )  ";
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
      }
   }

}
}
#endif
