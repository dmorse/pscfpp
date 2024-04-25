#ifndef RPC_LINEAR_RESPONSE_ANALYZER_TPP
#define RPC_LINEAR_RESPONSE_ANALYZER_TPP

#include "LinearResponseAnalyzer.h"

#include <rpc/simulate/Simulator.h>
#include <rpc/System.h>
#include <rpc/intra/IntraCorrelation.h>

#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>

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

//#include <unordered_map>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   LinearResponseAnalyzer<D>::LinearResponseAnalyzer(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      isInitialized_(false),
      kSize_(1)
   {  setClassName("LinearResponseAnalyzer"); }


   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void LinearResponseAnalyzer<D>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read(in,"nSamplePerBlock", nSamplePerBlock_);
      read(in,"StepSize", stepSize_);
   }
   
   /*
   * LinearResponseAnalyzer setup
   */
   template <int D>
   void LinearResponseAnalyzer<D>::setup() 
   {
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      
      //Allocate variables
      
      IntVec<D> const & dimensions = system().mesh().dimensions();
      
      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
            kSize_ *= dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
            kSize_ *= (dimensions[i]/2 + 1);
         }
      }
      
      w1_.allocate(nMonomer);
      w2_.allocate(nMonomer);
      
      for (int i = 0; i < nMonomer; ++i) {
         w1_[i].allocate(dimensions);
         w2_[i].allocate(dimensions);
      }
      
      intraCorrelation_.allocate(kMeshDimensions_);
      realError_.allocate(dimensions);
      realErrorDft_.allocate(dimensions);
      estimateDft_.allocate(dimensions);
      
      // Compute intraCorrelation (homopolymer)
      intraCorrelation_ = system().intraCorrelation().computeIntraCorrelations();
      
      nWave_ = kSize_;
      accumulators_.allocate(nWave_);
      
      // Clear accumulators
      for (int i = 0; i < nWave_; ++i) {
         accumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
         accumulators_[i].clear();
      }
      isInitialized_ = true;
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      
   }

   /* 
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void LinearResponseAnalyzer<D>::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         updateAccumulators();
      }
   }
   
   /*
   * Update accumulators for all current values.
   */
   template<int D>
   void LinearResponseAnalyzer<D>::updateAccumulators(){
      IntVec<D> const & dimensions = system().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j<  meshSize; j++){
            w1_[i][j] = system().w().rgrid(i)[j];
         }
      }
      
      for (int k = 0; k < meshSize; k++){
         //Random number generator
         double r = simulator().random().uniform(-stepSize_,stepSize_);
         for (int i = 0; i < nMonomer; i++){
            w2_[i][k] = w1_[i][k] + r;
         }
      }
      
      system().setWRGrid(w2_);
      
      // Compute real incompressibility error
      system().compute();
      for (int i = 0; i<  meshSize; i++){
         double real = 1;
         for (int j = 0; j <nMonomer ; j++){
            real -=  system().c().rgrid(j)[i];
         }
         realError_[i] = real;
      }
      
      system().fft().forwardTransform(realError_, realErrorDft_);
      
      RField<D> resid_;
      RFieldDft<D> residK_;
      const double vMonomer = system().mixture().vMonomer();

      resid_.allocate(dimensions);
      residK_.allocate(dimensions);

      // dw
      for (int i = 0; i < meshSize; ++i) {
         resid_[i] = w2_[0][i] - w1_[0][i];
      }
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      
      // Residual combine with Linear response factor (dphi = -vR* dw)
      MeshIterator<D> itr;
      itr.setDimensions(residK_.dftDimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         residK_[itr.rank()][0] *= vMonomer * intraCorrelation_[itr.rank()];
         residK_[itr.rank()][1] *= vMonomer * intraCorrelation_[itr.rank()];
      }
      
      // Convert back to real Space
      // system().fft().inverseTransform(residK_, resid_);
      
      // Define iterator
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      double sum = 0;
      
      // Compute sum of estimation over all fourier mode
      for (iter.begin(); !iter.atEnd(); ++iter) {
         std::complex<double> real(realErrorDft_[iter.rank()][0], realErrorDft_[iter.rank()][1]);
         std::complex<double> estimate(residK_[iter.rank()][0], residK_[iter.rank()][1]);
         sum += (estimate.real()* estimate.real() + estimate.imag()* estimate.imag())/kSize_;
      }
      
      for (iter.begin(); !iter.atEnd(); ++iter) {
         std::complex<double> real(realErrorDft_[iter.rank()][0], realErrorDft_[iter.rank()][1]);
         std::complex<double> estimate(residK_[iter.rank()][0], residK_[iter.rank()][1]);
         accumulators_[iter.rank()].sample(abs(real-estimate)/(std::sqrt(sum)));
         //accumulators_[iter.rank()].sample(abs((real-estimate)));
      }
      
      // Set back to original field
      system().setWRGrid(w1_);
      system().compute();
   
   }
   
  
   /*
   * Output final results to output file.
   */
   template <int D>  
   void LinearResponseAnalyzer<D>::output() 
   {
      // Output structure factors to one file
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "\t" << "k";
      outputFile_ << "\t" <<"\t" <<"phi_real/phi_LR";
      outputFile_ << "\t" <<"\t" <<"std";
      outputFile_<< std::endl;
      // Define iterator
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
         Gsq = system().unitCell().ksq(Gmin);
         
         outputFile_ << Dbl(sqrt(Gsq), 18,8) ;
         outputFile_ << Dbl(accumulators_[iter.rank()].average(),18,8);
         outputFile_ << Dbl(accumulators_[iter.rank()].stdDeviation(),18,8)<< std::endl;
      }
      outputFile_.close();
   }


}
}
#endif 
