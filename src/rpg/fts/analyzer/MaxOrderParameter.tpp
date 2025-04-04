#ifndef RPG_MAX_ORDER_PARAMETER_TPP
#define RPG_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/resources.h>

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <fftw3.h>

#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>

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
      UTIL_CHECK(system().w().hasData());
     
      // Precondition: Require that the system has two monomer types
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      AverageAnalyzer<D>::setup();

      if (!simulator().hasWc()){
         simulator().computeWc();
      }
   
      //Allocate variables
      IntVec<D> const & dimensions = system().mesh().dimensions();
      if (!isInitialized_){
         wc0_.allocate(dimensions);
         wK_.allocate(dimensions);
      }
      
      // Compute Fourier space dimension
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
            kSize_ *= dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
            kSize_ *= (dimensions[i]/2 + 1);
         }
      }
      
      isInitialized_ = true;
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }
   
   template <int D>
   double MaxOrderParameter<D>::compute() 
   {
      UTIL_CHECK(system().w().hasData());
   
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
