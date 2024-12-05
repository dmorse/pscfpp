#ifndef RPG_MAX_ORDER_PARAMETER_TPP
#define RPG_MAX_ORDER_PARAMETER_TPP

#include "MaxOrderParameter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>

#include <prdc/cuda/RField.h>

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
   MaxOrderParameter<D>::MaxOrderParameter(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      kSize_(1),
      hasAverage_(true),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("MaxOrderParameter"); }


   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void MaxOrderParameter<D>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOptional(in, "hasAverage", hasAverage_);
      readOutputFileName(in);
      readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
      
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "    chi       " << "MaxOrderParameter" << "\n";
   }
   
   /*
   * MaxOrderParameter setup
   */
   template <int D>
   void MaxOrderParameter<D>::setup() 
   {
      UTIL_CHECK(system().w().hasData());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      if (nMonomer != 2) {
         UTIL_THROW("The MaxOrderParameter Analyzer is designed specifically for diblock copolymer system. Please verify the number of monomer types in your system.");
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
      
      // Clear accumulators
      if (hasAverage_){
         accumulator_.clear();
      }
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }

   /* 
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void MaxOrderParameter<D>::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;
      computeMaxOrderParameter();
      
      if (hasAverage_){
         accumulator_.sample(maxOrderParameter_);
      }
      
      double chi =  system().interaction().chi(0,1);
      UTIL_CHECK(outputFile_.is_open());
      outputFile_ << Dbl(chi);
      outputFile_ << Dbl(maxOrderParameter_);
      outputFile_<< "\n";
   }
   
   template <int D>
   void MaxOrderParameter<D>::computeMaxOrderParameter()
   {
      const int meshSize = system().domain().mesh().size();
      RField<D> psi;
      psi.allocate(kMeshDimensions_);
      
      // GPU resources with meshSize threads
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Conver W_(r) to fourier mode W_(k)
      assignReal<<<nBlocks, nThreads>>>
            (wc0_.cArray(), simulator().wc(0).cArray(), meshSize);
      system().fft().forwardTransform(wc0_, wK_);
      
      // GPU resources with kSize threads
      ThreadGrid::setThreadsLogical(kSize_, nBlocks, nThreads);
      
      // Comput W_(k)^2
      squaredMagnitudeComplex<<<nBlocks,nThreads>>>
            (wK_.cArray(), psi.cArray(), kSize_);
            
      // Obtain max[W_(k)^2]
      maxOrderParameter_ = (double)gpuMaxAbs(psi.cArray(), kSize_);
   }
   
   /*
   * Output final results to output file.
   */
   template <int D>  
   void MaxOrderParameter<D>::output() 
   {
      if (hasAverage_){
         Log::file() << std::endl;
         Log::file() << "At chi = " << system().interaction().chi(0,1) << "\n";
         Log::file() << "Time average of the MaxOrderParameter is: "
                     << Dbl(accumulator_.average())
                     << " +- " << Dbl(accumulator_.blockingError(), 9, 2) 
                     << "\n";
      } 
      
   }

}
}
#endif 
