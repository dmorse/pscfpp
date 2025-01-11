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
      readOutputFileName(in);
      readOptional(in, "hasAverage", hasAverage_);
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
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      if (nMonomer != 2) {
         UTIL_THROW("The MaxOrderParameter Analyzer is designed specifically for diblock copolymer system. Please verify the number of monomer types in your system.");
      }
      
      //Allocate variables
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
      
      // Clear accumulators
      if (hasAverage_){
         accumulator_.clear();
      }
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }

   /* 
   * Compute a sampled max order parameter and update the accumulator.
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
      UTIL_CHECK(system().w().hasData());
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
