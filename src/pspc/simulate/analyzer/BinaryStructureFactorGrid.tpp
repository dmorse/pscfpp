#ifndef PSPC_BINARY_STRUCTURE_FACTOR_GRID_TPP
#define PSPC_BINARY_STRUCTURE_FACTOR_GRID_TPP

#include "BinaryStructureFactorGrid.h"
#include <pspc/simulate/McSimulator.h>
#include <pspc/System.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/format/Dbl.h>
#include <util/accumulators/Average.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <fftw3.h>

namespace Pscf {
namespace Pspc 
{
   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactorGrid<D>::BinaryStructureFactorGrid(McSimulator<D>& mcSimulator, System<D>& system) 
    : Analyzer<D>(),
      mcSimulatorPtr_(&mcSimulator),
      systemPtr_(&(mcSimulator.system())),
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
         wKGrid_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wKGrid_[i].allocate(dimensions);
         }
      }
      
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().forwardTransform(system().w().rgrid()[i], wKGrid_[i]);
      }
      
      nWave_ = wKGrid_[0].capacity();
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
      for (int i = 0; i < wm.capacity(); ++i) {
          wm[i] = (system().w().rgrid(0)[i] - system().w().rgrid(1)[i])/2;
      }
      RFieldDft<D> wk;
      wk.allocate(dimensions);
      system().fft().forwardTransform(wm, wk);
      // Convert real grid to KGrid format
      for (int k=0; k< wk.capacity(); k++) {
         std::complex<double> wmKGrid(wk[k][0], wk[k][1]);
         double squared_magnitude = std::norm(wmKGrid);
         accumulators_[k].sample(squared_magnitude);
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
      itr.setDimensions(wKGrid_[0].dftDimensions());
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
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().forwardTransform(system().w().rgrid()[i], wKGrid_[i]);
      }
      IntVec<D> G; IntVec<D> Gmin; 
      MeshIterator<D> itr(wKGrid_[0].dftDimensions());
      double kSq;
      std::vector<double> qRoList(wKGrid_[0].capacity());
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
      outputFile_ << "\t" << "qRO";
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
