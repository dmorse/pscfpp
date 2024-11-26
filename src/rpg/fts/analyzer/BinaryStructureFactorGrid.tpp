#ifndef RPG_BINARY_STRUCTURE_FACTOR_GRID_TPP
#define RPG_BINARY_STRUCTURE_FACTOR_GRID_TPP

#include "BinaryStructureFactorGrid.h"

#include <prdc/cuda/HostField.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>
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

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactorGrid<D>::BinaryStructureFactorGrid(Simulator<D>& simulator, System<D>& system) 
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      isInitialized_(false)
   {  setClassName("BinaryStructureFactorGrid"); }

   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void BinaryStructureFactorGrid<D>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read(in,"nSamplePerBlock", nSamplePerBlock_);
   }
   
   /*
   * BinaryStructureFactorGrid setup
   */
   template <int D>
   void BinaryStructureFactorGrid<D>::setup() 
   {
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      if (nMonomer != 2) {
         UTIL_THROW("The BinaryStructureFactorGrid Analyzer is designed specifically for diblock copolymer system. Please verify the number of monomer types in your system.");
      }
      
      //Allocate variables
      IntVec<D> const & dimensions = system().mesh().dimensions();
      if (!isInitialized_){
         wkGrid_.allocate(nMonomer);
         wrGrid_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wkGrid_[i].allocate(dimensions);
            wrGrid_[i].allocate(dimensions);
         }
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
         updateAccumulators();
      }
   }
   
   /*
   * Update accumulators for all current values.
   */
   template <int D>   
   void BinaryStructureFactorGrid<D>::updateAccumulators() 
   {
      IntVec<D> const & dimensions = system().mesh().dimensions();
      RField<D> wm;
      wm.allocate(dimensions);
      
      // Compute W-
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(wm.capacity(), nBlocks, nThreads);
      
      pointWiseBinarySubtract<<<nBlocks, nThreads>>>
            (system().w().rgrid(0).cField(), system().w().rgrid(1).cField(), wm.cField(), wm.capacity());
      scaleReal<<<nBlocks, nThreads>>>(wm.cField(), 0.5, wm.capacity());
      
      // Convert real grid to KGrid format
      RFieldDft<D> wk;
      wk.allocate(dimensions);
      system().fft().forwardTransform(wm, wk);
      
      std::vector<std::complex<double>> wkCpu(kSize_);
      cudaMemcpy(wkCpu.data(), wk.cField(), kSize_ * sizeof(cudaComplex), cudaMemcpyDeviceToHost);
      for (int k=0; k< wk.capacity(); k++) {
         accumulators_[k].sample(norm(wkCpu[k]));
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
         structureFactors_[itr.rank()] = n / (chi * chi) * accumulators_[itr.rank()].average() - 1.0/(2.0*chi);
      }
   }
   
   // Average S(k) over k of equal magnitude
   template <int D>   
   void BinaryStructureFactorGrid<D>::averageStructureFactor() 
   {
      // Convert real grid to KGrid format
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> G; IntVec<D> Gmin; 
      MeshIterator<D> itr(kMeshDimensions_);
      double kSq;
      std::vector<double> qRoList(kSize_);
      computeRoSquare();
      for (itr.begin(); !itr.atEnd(); ++itr){
         // Obtain square magnitude of reciprocal basis vector
         G = itr.position();
         Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
         kSq = system().unitCell().ksq(Gmin);
         qRoList[itr.rank()] = sqrt(kSq * roSquare_);
      }
      
      UTIL_CHECK(qRoList.capacity() == structureFactors_.capacity());
      std::map<double, double> SMap;
      for (int i = 0; i < qRoList.capacity(); ++i) {
        double qRo = qRoList[i];
        double  s = structureFactors_[i];
        SMap[qRo] += s;
      }
      
      // Average structure factor with same magnitude value of qRo
      for (auto& i : SMap) {
        double qRo = i.first;
        double sum = i.second;
        int count = std::count(qRoList.begin(), qRoList.end(), qRo);
        i.second = sum / count;
      }
      
      // Average structure factor for qRo in range of +-tolerance
      double tolerance  = 1e-5;
      auto it = SMap.begin();
      double currentqRo = it->first;
      double sumS = it->second;
      int count = 1;
      for (++it; it != SMap.end(); ++it) {
         double key = it->first;
         double s = it->second;
         if (std::abs(key - currentqRo) <= tolerance) {
            sumS += s;
            count++;
         } else {
            averageSMap_[currentqRo] = sumS / count;
            currentqRo = key;
            sumS = s;
            count = 1;
         }
      }
      // Add last batch of data (last number within tolerance)
      averageSMap_[currentqRo] = sumS / count;
   }
   
   /**
   * Compute radius of gyration Ro^2 = Nb^2
   */
   template <int D>
   void BinaryStructureFactorGrid<D>::computeRoSquare()
   {
      int np = system().mixture().nPolymer();
      if (np > 0) {
         double nb; // number of monomer in each polymer
         double length; // fraction of monomer
         double kuhn; // statistical segment length of monomer
         int i; // molecule index
         int j; // block index
         for (i = 0; i < np; ++i) {
            nb = system().mixture().polymer(i).nBlock();
            roSquare_ = 0;
            for (j = 0; j < nb; ++j) {
               Block<D>& block = system().mixture().polymer(i).block(j);
               kuhn = block.kuhn();
               length = block.length();
               roSquare_ += length * kuhn * kuhn;
            }
         }
      }      
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
      outputFile_ << "\t" << "kRO";
      outputFile_ << "\t" <<"\t" <<"S(k)/(\u03C1 N)";
      outputFile_<< std::endl;
      for (const auto& i : averageSMap_) {
         double qRo = i.first;
         double averageS = i.second;
         outputFile_ << Dbl(qRo, 18, 8);
         outputFile_ << Dbl(averageS, 18, 8);
         outputFile_ << std::endl;
      }
      outputFile_.close();
   }

}
}
#endif 
