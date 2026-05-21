#ifndef RPC_BINARY_STRUCTURE_FACTOR_TPP
#define RPC_BINARY_STRUCTURE_FACTOR_TPP

#include "BinaryStructureFactor.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/FFT.h>
#include <prdc/cpu/WaveList.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/complex.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/global.h>

#include <iostream>
#include <complex>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactor<D>::BinaryStructureFactor(
                                       Simulator<D>& simulator,
                                       System<D>& system)
    : Analyzer<D>(simulator, system),
      kMeshDimensions_(),
      nWave_(0),
      nBunch_(0),
      keepWaveData_(false),
      isInitialized_(false)
   {
      ParamComposite::setClassName("BinaryStructureFactor");
      AnalyzerT::setFileMaster(system.fileMaster());
   }

   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void BinaryStructureFactor<D>::readParameters(std::istream& in)
   {
      // Precondition: Require that the system has two monomer types
      UTIL_CHECK(system().mixture().nMonomer() == 2);

      AnalyzerT::readInterval(in);
      AnalyzerT::readOutputFileName(in);
      keepWaveData_ = false;
      ParamComposite::readOptional<bool>(in, "keepWaveData", 
                                         keepWaveData_);
      isInitialized_ = true;
   }

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D>::setup()
   {
      UTIL_CHECK(isInitialized_);

      // Compute and compute mesh dimensions
      Mesh<D> const & mesh = system().domain().mesh();
      IntVec<D> const & rMeshDimensions = mesh.dimensions();
      FFT<D>::computeKMesh(rMeshDimensions, kMeshDimensions_, nWave_);

      // If needed, allocate arrays indexed by wave id
      if (!wm_.isAllocated()){
         UTIL_CHECK(!wk_.isAllocated());
         UTIL_CHECK(!waveBunchIds_.isAllocated());
         UTIL_CHECK(!waveWeights_.isAllocated());
         wm_.allocate(rMeshDimensions);
         wk_.allocate(rMeshDimensions);
         waveBunchIds_.allocate(nWave_);
         waveWeights_.allocate(nWave_);
         if (keepWaveData_) {
            UTIL_CHECK(!waveAccumulators_.isAllocated());
            waveAccumulators_.allocate(nWave_);
         }
      }
      UTIL_CHECK(wm_.capacity() == mesh.size());
      UTIL_CHECK(wk_.capacity() == nWave_);
      UTIL_CHECK(waveBunchIds_.capacity() == nWave_);
      UTIL_CHECK(waveWeights_.capacity() == nWave_);
      if (keepWaveData_) {
         UTIL_CHECK(waveAccumulators_.capacity() == nWave_);
      }

      // Sort waves and set nBunch
      WaveList<D>& waveList = system().waveList();
      if (!waveList.hasKSq()) {
         waveList.computeKSq();
      }
      waveList.sortWaves();
      nBunch_ = waveList.nBunch();
      UTIL_CHECK(nBunch_ > 0);

      // If needed, allocate or re-allocate arrays indexed by bunch id
      if (bunchAccumulators_.isAllocated()) {
         int const m = bunchAccumulators_.capacity();
         UTIL_CHECK(bunchWavenumbers_.capacity() == m);
         UTIL_CHECK(bunchValues_.capacity() == m);
         if (bunchAccumulators_.capacity() < nBunch_) {
            bunchAccumulators_.deallocate();
            bunchWavenumbers_.deallocate();
            bunchValues_.deallocate();
         }
      }
      if (!bunchAccumulators_.isAllocated()) {
         UTIL_CHECK(!bunchWavenumbers_.isAllocated());
         UTIL_CHECK(!bunchValues_.isAllocated());
         bunchAccumulators_.allocate(nBunch_);
         bunchWavenumbers_.allocate(nBunch_);
         bunchValues_.allocate(nBunch_);
      }
      int m = bunchAccumulators_.capacity();
      UTIL_CHECK(m >= nBunch_);
      UTIL_CHECK(bunchWavenumbers_.capacity() == m);
      UTIL_CHECK(bunchValues_.capacity() == m);

      // Initialize empty arrays
      for (int ib = 0; ib < m; ++ib) {
         bunchWavenumbers_[ib] = 0.0;
         bunchValues_[ib] = 0.0;
         bunchAccumulators_[ib].clear();
      }
      for (int iw = 0; iw < nWave_; ++iw) {
         waveBunchIds_[iw] = -1;
         waveWeights_[iw] = 0.0;
      }
      if (keepWaveData_) {
         for (int iw = 0; iw < nWave_; ++iw) {
            waveAccumulators_[iw].clear();
         }
      }

      // Define references to WaveList data structures
      Array< double > const & kSq = waveList.kSq();
      Array< int > const & sortedIds = waveList.sortedIds();
      Array<bool> const & implicit = waveList.implicitInverse();
      GArray< Pair<int> > const & bunches = waveList.sortedBunches();
      UTIL_CHECK(kSq.capacity() == nWave_);
      UTIL_CHECK(sortedIds.capacity() == nWave_);
      UTIL_CHECK(implicit.capacity() == nWave_);
      UTIL_CHECK(bunches.size() == nBunch_);

      // Set bunchWavenumbers_, waveBunchIds_, and waveWeights_
      double wavenumberSq, newWavenumber, oldWavenumber;
      double count, sum, tot;
      int begin, end, k, iw;
      int nw = 0;
      for (int ib = 0; ib < nBunch_; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];

         // Compute bunchWavenumbers_[ib] from first wave in bunch
         iw = sortedIds[begin]; // wave id of first wave in bunch
         wavenumberSq = kSq[iw];
         UTIL_CHECK(wavenumberSq >= 0.0);
         newWavenumber = std::sqrt(wavenumberSq);
         bunchWavenumbers_[ib] = newWavenumber;
         if (ib > 0) {
            UTIL_CHECK(newWavenumber - oldWavenumber > 1.0E-8);
         }
         oldWavenumber = newWavenumber;

         // Set bunchIds, unnormalized weights for all waves in this bunch
         sum = 0.0;
         for (k = begin; k < end; ++k) {
            iw = sortedIds[k];
            // Check that each wave is treated only once
            UTIL_CHECK(waveBunchIds_[iw] == -1);
            UTIL_CHECK(std::abs(waveWeights_[iw]) < 1.0E-8);
            waveBunchIds_[iw] = ib;
            if (implicit[iw]) {
               count = 2.0;
            } else {
               count = 1.0;
            }
            waveWeights_[iw] = count;
            sum += count;
            UTIL_CHECK(std::abs(kSq[iw] - wavenumberSq) < 1.0E-8);
         }

         // Normalize weights for all waves in this bunch
         tot = 0.0;
         for (k = begin; k < end; ++k) {
            iw = sortedIds[k];
            waveWeights_[iw] /= sum;
            tot += waveWeights_[iw];
            ++nw;
         }
         UTIL_CHECK(std::abs(tot - 1.0) < 1.0E-8);

      }
      UTIL_CHECK(nw == nWave_);

   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D>::sample(long iStep)
   {
      if (AnalyzerT::isAtInterval(iStep)) {

         // Preconditions
         UTIL_CHECK(isInitialized_);
         UTIL_CHECK(nWave_ > 0);
         UTIL_CHECK(nBunch_ > 0);
         UTIL_CHECK(nWave_ >= nBunch_);
         UTIL_CHECK(wk_.capacity() == nWave_);
         UTIL_CHECK(waveBunchIds_.capacity() == nWave_);
         UTIL_CHECK(waveWeights_.capacity() == nWave_);
         int m = bunchAccumulators_.capacity();
         UTIL_CHECK(m >= nBunch_);
         UTIL_CHECK(bunchWavenumbers_.capacity() >= m);
         UTIL_CHECK(bunchValues_.capacity() >= m);
         UTIL_CHECK(system().w().hasData());

         // Compute W_{-}(r)
         RField<D> const & wa = system().w().rgrid(0);
         RField<D> const & wb = system().w().rgrid(1);
         VecOp::subVV(wm_, wa, wb);
         VecOp::mulEqS(wm_, 0.5);

         // Fourier transform W_{-}(r)
         system().domain().fft().forwardTransform(wm_, wk_);

         // Initialize bunch average values to zero
         for (int ib = 0; ib < nBunch_; ++ib) {
            bunchValues_[ib] = 0.0;
         }

         // Set constant coefficients
         double const vSystem  = system().domain().unitCell().volume();
         double const vMonomer = system().mixture().vMonomer();
         double const chi = system().interaction().chi(0,1);
         double a_ = vSystem / (chi * chi * vMonomer * vMonomer);
         double b_ = 0.5 / (chi * vMonomer);

         // Compute structure factors for all waves, add to bunchValues_ 
         double value;
         int ib;
         for (int iw = 0; iw < nWave_; ++iw) {
            value = a_ * absSq( wk_[iw] );
            value -= b_;
            if (keepWaveData_) {
               waveAccumulators_[iw].sample(value);
            }
            ib = waveBunchIds_[iw];
            bunchValues_[ib] += waveWeights_[iw] * value;
         }

         // Pass spherically averaged values to bunchAccumulators_
         for (ib = 0; ib < nBunch_; ++ib) {
            bunchAccumulators_[ib].sample(bunchValues_[ib]);
         }

      }
   }

   /*
   * Output final results to output file.
   */
   template <int D>
   void BinaryStructureFactor<D>::output()
   {
      std::string outputFileName;
      std::ofstream outputFile_;

      // Output spherical average values of S(q) for bunches
      outputFileName = AnalyzerT::outputFileName();
      if (keepWaveData_) {
         outputFileName += "_ave";
      }
      AnalyzerT::fileMaster().openOutputFile(outputFileName, outputFile_);
      for (int i = 0; i < nBunch_; ++i) {
         outputFile_ << Dbl(bunchWavenumbers_[i], 18, 8);
         outputFile_ << Dbl(bunchAccumulators_[i].average(), 18, 8);
         outputFile_ << std::endl;
      }
      outputFile_.close();

      // Optionally output S(q) values for individual waves
      if (keepWaveData_) {
         outputFileName = AnalyzerT::outputFileName();
         outputFileName += "_wave";
         AnalyzerT::fileMaster().openOutputFile(outputFileName, 
                                                outputFile_);
         MeshIterator<D> iter(kMeshDimensions_);
         IntVec<D> p; 
         int iw, j;
         for (iter.begin(); !iter.atEnd(); ++iter) {
            iw = iter.rank();
            p = iter.position();
            outputFile_ << Int(iw)
                        << Dbl(waveAccumulators_[iw].average(), 18, 8);
            for (j = 0; j < D; ++j) {
               outputFile_ << Int(p[j], 5);
            }
            outputFile_ << std::endl;
         }
         outputFile_.close();
      }

   }

}
}
#endif
