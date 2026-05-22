#ifndef RP_BINARY_STRUCTURE_FACTOR_TPP
#define RP_BINARY_STRUCTURE_FACTOR_TPP

#include "BinaryStructureFactor.h"

#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/containers/Array.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/global.h>

#include <iostream>
#include <complex>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   BinaryStructureFactor<D,T>::BinaryStructureFactor(
                                typename T::Simulator& simulator,
                                typename T::System& system)
    : AnalyzerT(simulator, system),
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
   template <int D, class T>
   void BinaryStructureFactor<D,T>::readParameters(std::istream& in)
   {
      // Precondition: Require that the system has two monomer types
      UTIL_CHECK(system().mixture().nMonomer() == 2);

      AnalyzerT::readInterval(in);
      AnalyzerT::readOutputFileName(in);
      keepWaveData_ = false;
      ParamComposite::readOptional(in, "keepWaveData", keepWaveData_);
      isInitialized_ = true;
   }

   /*
   * Allocate memory arrays with dimensions that depend only on mesh.
   */
   template <int D, class T>
   void BinaryStructureFactor<D,T>::allocate()
   {
      UTIL_CHECK(isInitialized_);

      // Compute and compute mesh dimensions
      Mesh<D> const & mesh = system().domain().mesh();
      IntVec<D> const & rMeshDimensions = mesh.dimensions();
      FFTT::computeKMesh(rMeshDimensions, kMeshDimensions_, nWave_);

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
   }

   /*
   * Allocate and initialize data structures that involve wave bunches.
   */
   template <int D, class T>
   void BinaryStructureFactor<D,T>::findWaveBunches(
                                  Array<double> const & kSq,
                                  Array<bool> const & implicit)
   {
      UTIL_CHECK(kSq.capacity() == nWave_);
      UTIL_CHECK(implicit.capacity() == nWave_);

      // Sort waves and set nBunch
      typename T::WaveList& waveList = system().waveList();
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
      Array< int > const & sortedIds = waveList.sortedIds();
      GArray< Pair<int> > const & bunches = waveList.sortedBunches();
      UTIL_CHECK(sortedIds.capacity() == nWave_);
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
            // Check that each wave is only processed once
            UTIL_CHECK(waveBunchIds_[iw] == -1);
            UTIL_CHECK(std::abs(waveWeights_[iw]) < 1.0E-8);
            waveBunchIds_[iw] = ib;
            if (implicit[iw]) {
               count = 2.0;
            } else {
               count = 1.0;
            }
            // count = 1.0;   // Uncomment to weight waves equally
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
   * Compute W_{-} and it Fourier transform.
   */
   template <int D, class T>
   void BinaryStructureFactor<D,T>::computeW()
   {
      // Preconditions
      UTIL_CHECK(isInitialized_);
      UTIL_CHECK(nWave_ > 0);
      UTIL_CHECK(wk_.capacity() == nWave_);

      // Compute W_{-}(r)
      typename T::RField const & wa = system().w().rgrid(0);
      typename T::RField const & wb = system().w().rgrid(1);
      VecOp::subVV(wm_, wa, wb);
      VecOp::mulEqS(wm_, 0.5);

      // Fourier transform W_{-}(r)
      system().domain().fft().forwardTransform(wm_, wk_);
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D, class T>
   void BinaryStructureFactor<D,T>::computeS(
                                    Array<typename T::Complex> const & wk)
   {
      // Preconditions
      UTIL_CHECK(isInitialized_);
      UTIL_CHECK(nBunch_ > 0);
      UTIL_CHECK(nWave_ >= nBunch_);
      int m = bunchAccumulators_.capacity();
      UTIL_CHECK(m >= nBunch_);
      UTIL_CHECK(bunchWavenumbers_.capacity() >= m);
      UTIL_CHECK(bunchValues_.capacity() >= m);
      UTIL_CHECK(waveBunchIds_.capacity() == nWave_);
      UTIL_CHECK(waveWeights_.capacity() == nWave_);
      UTIL_CHECK(wk.capacity() == nWave_);

      // Initialize bunch average values to zero
      for (int ib = 0; ib < nBunch_; ++ib) {
         bunchValues_[ib] = 0.0;
      }

      // Set coefficients a_ and b_
      double const vSystem  = system().domain().unitCell().volume();
      double const vMonomer = system().mixture().vMonomer();
      double const chi = system().interaction().chi(0,1);
      double a_ = vSystem / (chi * chi * vMonomer * vMonomer);
      double b_ = 0.5 / (chi * vMonomer);

      // Compute structure factors for all waves, add to bunchValues_
      double value;
      int ib;
      for (int iw = 0; iw < nWave_; ++iw) {
         value = a_ * absSq( wk[iw] );
         value -= b_;
         if (keepWaveData_) {
            waveAccumulators_[iw].sample(value);
         }
         ib = waveBunchIds_[iw];
         bunchValues_[ib] += waveWeights_[iw] * value;
      }

      // Pass bunchValues_ to bunchAccumulators_
      for (ib = 0; ib < nBunch_; ++ib) {
         bunchAccumulators_[ib].sample(bunchValues_[ib]);
      }

   }

   /*
   * Output final results to output file.
   */
   template <int D, class T>
   void BinaryStructureFactor<D,T>::output()
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
