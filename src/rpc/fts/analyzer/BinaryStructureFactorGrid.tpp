#ifndef RPC_BINARY_STRUCTURE_FACTOR_GRID_TPP
#define RPC_BINARY_STRUCTURE_FACTOR_GRID_TPP

#include "BinaryStructureFactorGrid.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>
#include <complex>
#include <vector>
#include <map>
#include <algorithm>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactorGrid<D>::BinaryStructureFactorGrid(
                                          Simulator<D>& simulator,
                                          System<D>& system)
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      isInitialized_(false),
      nSamplePerBlock_(1),
      kMeshDimensions_(0),
      kSize_(0)
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
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Compute Fourier space kMeshDimensions_ and kSize_
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      if (!isInitialized_){
         wKGrid_.allocate(nMonomer);
         wm_.allocate(dimensions);
         wk_.allocate(dimensions);
         for (int i = 0; i < nMonomer; ++i) {
            wKGrid_[i].allocate(dimensions);
         }

         // Compute qList
         qList_.resize(kSize_);
         IntVec<D> G;
         IntVec<D> Gmin;
         UnitCell<D> const & unitCell = system().domain().unitCell();
         MeshIterator<D> itr(kMeshDimensions_);
         double Gsq;
         for (itr.begin(); !itr.atEnd(); ++itr){
            // Obtain square magnitude of reciprocal basis vector
            G = itr.position();
            Gmin = shiftToMinimum(G, dimensions, unitCell);
            Gsq = unitCell.ksq(Gmin);
            qList_[itr.rank()] = sqrt(Gsq);
         }

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

         // Compute W-
         for (int i = 0; i < wm_.capacity(); ++i) {
             wm_[i] = (system().w().rgrid(0)[i] - system().w().rgrid(1)[i])/2;
         }

         // Convert real grid to KGrid format
         system().domain().fft().forwardTransform(wm_, wk_);

         // Pass square magnitudes of Fourier components to accumulators
         for (int k=0; k< wk_.capacity(); k++) {
            std::complex<double> wmKGrid(wk_[k][0], wk_[k][1]);
            double squared_magnitude = std::norm(wmKGrid);
            accumulators_[k].sample(squared_magnitude);
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
      itr.setDimensions(wKGrid_[0].dftDimensions());
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
      // Convert real grid to KGrid format
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         system().domain().fft().forwardTransform(system().w().rgrid()[i],
                                                  wKGrid_[i]);
      }

      std::map<double, double> SMap;
      {
         int qListcapacity = (int)qList_.capacity();
         double q, s;
         UTIL_CHECK(qListcapacity == structureFactors_.capacity());
         for (int i = 0; i < qListcapacity; ++i) {
           q = qList_[i];
           s = structureFactors_[i];
           SMap[q] += s;
         }
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
