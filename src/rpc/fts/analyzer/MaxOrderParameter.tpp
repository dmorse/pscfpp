#ifndef RPC_MAX_ORDER_PARAMETER_TPP
#define RPC_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/System.h>

#include <prdc/cpu/RField.h>

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
      AverageAnalyzer<D>::setup();
      
      // Allocate variables
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      if (!isInitialized_){
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
  
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);
      
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      
      MeshIterator<D> itr;
      itr.setDimensions(kMeshDimensions_);
      std::vector<double> psi(kSize_);
      
      // Conver W_(r) to fourier mode W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);
      
      for (itr.begin(); !itr.atEnd(); ++itr) {
         std::complex<double> wK(wK_[itr.rank()][0], wK_[itr.rank()][1]);
         psi[itr.rank()] = std::norm(wK);
      }
      
      auto maxOrderParameterPtr = std::max_element(psi.begin(), psi.end());
      maxOrderParameter_ = *maxOrderParameterPtr;
      
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
