#ifndef RPG_MAX_ORDER_PARAMETER_TPP
#define RPG_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/resources.h>

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

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
      // Local copies of data
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().mesh().dimensions();

      // Preconditions
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(nMonomer == 2);  // Require system has 2 monomer types

      AverageAnalyzer<D>::setup();
      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      // Compute Fourier space dimension kMeshDimensions_ and kSize_
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wc0_.allocate(dimensions);
         wK_.allocate(dimensions);
      }

      isInitialized_ = true;
   }

   template <int D>
   double MaxOrderParameter<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());

      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      const int meshSize = system().domain().mesh().size();
      RField<D> psi;
      psi.allocate(kMeshDimensions_);

      // Conver W_(r) to fourier mode W_(k)
      VecOp::eqV(wc0_, simulator().wc(0));
      system().fft().forwardTransform(wc0_, wK_);

      // Comput W_(k)^2
      VecOp::sqNormV(psi, wK_);

      // Obtain max[W_(k)^2]
      maxOrderParameter_ = Reduce::maxAbs(psi);

      return maxOrderParameter_;
   }

   template <int D>
   void MaxOrderParameter<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double chi= system().interaction().chi(0,1);

         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(chi);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }

}
}
#endif
