#ifndef RPG_BINARY_STRUCTURE_FACTOR_GRID_TPP
#define RPG_BINARY_STRUCTURE_FACTOR_GRID_TPP

#include "BinaryStructureFactorGrid.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/cuda/resources.h>
#include <prdc/cuda/complex.h>

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>

#include <fftw3.h>

#include <vector>
#include <unordered_map>
#include <algorithm>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactorGrid<D>::BinaryStructureFactorGrid(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      isInitialized_(false),
      nSamplePerBlock_(1)
   {  setClassName("BinaryStructureFactorGrid"); }

   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void BinaryStructureFactorGrid<D>::readParameters(std::istream& in) 
   {
      // Precondition: Require that the system has two monomer types
      UTIL_CHECK(system().mixture().nMonomer() == 2);

      readInterval(in);
      readOutputFileName(in);
      readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
   }
   
   /*
   * BinaryStructureFactorGrid setup
   */
   template <int D>
   void BinaryStructureFactorGrid<D>::setup() 
   {
      //Allocate variables
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      if (!isInitialized_){
         wm_.allocate(dimensions);
         wk_.allocate(dimensions);
      }
      
      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
         }
      }
      kSize_ = 1;
      for(int i = 0; i < D; ++i) {
         kSize_ *= kMeshDimensions_[i];
      }
      
      nWave_ = kSize_;
      structureFactors_.allocate(nWave_);
      accumulators_.allocate(nWave_);
      
      // Convert real grid to KGrid format
      qList_.resize(kSize_);
      RField<D> const & kSq = system().domain().waveList().kSq();
      HostDArray<cudaReal> kSqHost(kSize_);
      kSqHost = kSq;
      MeshIterator<D> itr;
      itr.setDimensions(kMeshDimensions_);
      for (itr.begin(); !itr.atEnd(); ++itr){
         qList_[itr.rank()] = sqrt(double(kSqHost[itr.rank()]));
      }
      
      isInitialized_ = true;
      
      // Clear accumulators
      for (int i = 0; i < nWave_; ++i) {
         structureFactors_[i] = 0.0;
         accumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
         accumulators_[i].clear();
      }
      
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
   }

   /* 
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void BinaryStructureFactorGrid<D>::sample(long iStep) 
   {
      UTIL_CHECK(system().w().hasData());
      
      if (isAtInterval(iStep))  {
         IntVec<D> const & dimensions = system().mesh().dimensions();
         
         // Compute W: (rgrid(0) - rgrid(1)) / 2
         VecOp::addVcVc(wm_, system().w().rgrid(0), 0.5, 
                        system().w().rgrid(1), -0.5);
         
         // Convert real grid to KGrid format
         system().fft().forwardTransform(wm_, wk_);
         
         HostDArray<cudaComplex> wkCpu(kSize_);
         wkCpu = wk_; // copy from device to host
         for (int k = 0; k < wk_.capacity(); k++) {
            accumulators_[k].sample(absSq<cudaComplex, cudaReal>(wkCpu[k]));
         }
      }
      
   }
   
   template <int D>
   void BinaryStructureFactorGrid<D>::computeStructureFactor()
   {
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      double n = vSystem / vMonomer;
      double chi= system().interaction().chi(0,1);
      MeshIterator<D> itr;
      itr.setDimensions(kMeshDimensions_);
      for (itr.begin(); !itr.atEnd(); ++itr) {
         // Compute vS(q)
         structureFactors_[itr.rank()] = n / (chi * chi) * accumulators_[itr.rank()].average() - 1.0/(2.0*chi);
         
         // Compute S(q)
         structureFactors_[itr.rank()] /= vMonomer;
      }
   }
   
   // Average S(k) over k of equal magnitude
   template <int D>   
   void BinaryStructureFactorGrid<D>::averageStructureFactor() 
   {
      UTIL_CHECK(qList_.capacity() == structureFactors_.capacity());
      
      std::map<double, double> SMap;
      for (int i = 0; i < qList_.capacity(); ++i) {
        double q = qList_[i];
        double  s = structureFactors_[i];
        SMap[q] += s;
      }
      
      // Average structure factor with same magnitude value of q
      for (auto& i : SMap) {
        double q = i.first;
        double sum = i.second;
        int count = std::count(qList_.begin(), qList_.end(), q);
        i.second = sum / count;
      }
      
      // Average structure factor for q in range of +-tolerance
      double tolerance  = 1e-5;
      auto it = SMap.begin();
      double currentq = it->first;
      double sumS = it->second;
      int count = 1;
      for (++it; it != SMap.end(); ++it) {
         double key = it->first;
         double s = it->second;
         if (std::abs(key - currentq) <= tolerance) {
            sumS += s;
            count++;
         } else {
            averageSMap_[currentq] = sumS / count;
            currentq = key;
            sumS = s;
            count = 1;
         }
      }
      
      // Add last batch of data (last number within tolerance)
      averageSMap_[currentq] = sumS / count;
   }
   

      
   
   /*
   * Output final results to output file.
   */
   template <int D>  
   void BinaryStructureFactorGrid<D>::output() 
   {
      computeStructureFactor();
      averageStructureFactor();

      // Output structure factors to one file
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "\t" << "q";
      outputFile_ << "\t" <<"\t" <<"S(q)/(\u03C1 N)";
      outputFile_<< std::endl;
      for (const auto& i : averageSMap_) {
         double q = i.first;
         double averageS = i.second;
         outputFile_ << Dbl(q, 18, 8);
         outputFile_ << Dbl(averageS, 18, 8);
         outputFile_ << std::endl;
      }
      outputFile_.close();
   }

}
}
#endif 
