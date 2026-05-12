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

#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/complex.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/global.h>

#include <iostream>

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
      kMeshDimensions_(0),
      kSize_(0),
      nBunch_(0),
      nSamplePerBlock_(1),
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
      ParamComposite::readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
      isInitialized_ = true;
   }

   /*
   * BinaryStructureFactor setup
   */
   template <int D>
   void BinaryStructureFactor<D>::setup()
   {
      UTIL_CHECK(isInitialized_);

      // Store and/or compute mesh dimensions
      Mesh<D> const & mesh = system().domain().mesh();
      IntVec<D> const & rMeshDimensions = mesh.dimensions();
      int const & rSize = mesh.size();
      FFT<D>::computeKMesh(rMeshDimensions, kMeshDimensions_, kSize_);

      // As needed, allocate arrays indexed by wave id
      if (!wm_.isAllocated()){
         UTIL_CHECK(!wk_.isAllocated());
         UTIL_CHECK(!bunchIds_.isAllocated());
         UTIL_CHECK(!weights_.isAllocated());
         wm_.allocate(rMeshDimensions);
         wk_.allocate(rMeshDimensions);
         bunchIds_.allocate(kSize_);
         weights_.allocate(kSize_);
      }
      UTIL_CHECK(wm_.capacity() == rSize);
      UTIL_CHECK(wk_.capacity() == kSize_);
      UTIL_CHECK(bunchIds_.capacity() == kSize_);
      UTIL_CHECK(weights_.capacity() == kSize_);

      // Sort waves and set nBunch
      WaveList<D>& waveList = system().waveList();
      if (!waveList.hasKSq()) {
         waveList.computeKSq();
      }
      waveList.sortWaves();
      nBunch_ = waveList.nBunch();

      // As needed, allocate arrays indexed by bunch id
      if (accumulators_.isAllocated()) {
         int const m = accumulators_.capacity();
         UTIL_CHECK(wavenumbers_.capacity() == m);
         UTIL_CHECK(values_.capacity() == m);
         if (accumulators_.capacity() < nBunch_) {
            accumulators_.deallocate();
            wavenumbers_.deallocate();
            values_.deallocate();
         }
      }
      if (!accumulators_.isAllocated()) {
         UTIL_CHECK(!wavenumbers_.isAllocated());
         UTIL_CHECK(!values_.isAllocated());
         accumulators_.allocate(nBunch_);
         wavenumbers_.allocate(nBunch_);
         values_.allocate(nBunch_);
      }
      int n = accumulators_.capacity();
      UTIL_CHECK(n >= nBunch_);
      UTIL_CHECK(wavenumbers_.capacity() == n);
      UTIL_CHECK(values_.capacity() == n);

      // Initialize empty arrays
      int ib;
      for (ib = 0; ib < n; ++ib) {
         accumulators_[ib].setNSamplePerBlock(nSamplePerBlock_);
         accumulators_[ib].clear();
         wavenumbers_[ib] = 0.0;
         values_[ib] = 0.0;
      }
      int iw;
      for (iw = 0; iw < kSize_; ++iw) {
         weights_[iw] = 0.0;
      }

      // Set values for wavenumbers_, bunchIds_, and weights_
      Array< double > const & kSq = waveList.kSq();
      Array< int > const & sortedIds = waveList.sortedIds();
      GArray< Pair<int> > const & bunches = waveList.sortedBunches();
      Array<bool> const & implicit = waveList.implicitInverse();
      double wavenumberSq, count, sum, tot;
      int begin, end, k;
      for (ib = 0; ib < nBunch_; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];
         iw = sortedIds[begin];

         // Set wavenumber for this bunch
         wavenumberSq = kSq[iw];
         UTIL_CHECK(wavenumberSq >= 0.0);
         wavenumbers_[ib] = std::sqrt(wavenumberSq);

         // Set bunchIds and unnormalized weights for waves in this bunch
         sum = 0.0;
         for (k = begin; k < end; ++k) {
            iw = sortedIds[k];
            bunchIds_[iw] = ib;
            if (implicit[iw]) {
               count = 2.0;
            } else {
               count = 1.0;
            }
            weights_[iw] = count;
            sum += count;
         }
         #if 0
         std::cout << std::endl << Int(ib)  << "  "
                   << Dbl(wavenumbers_[ib])
                   << Dbl(wavenumberSq) << Dbl(sum);
         #endif

         // Normalize weights for waves in this bunch
         tot = 0.0;
         for (k = begin; k < end; ++k) {
            iw = sortedIds[k];
            weights_[iw] /= sum;
            tot += weights_[iw];
            #if 0
            IntVec<D> const & minImage = waveList.minImages()[iw];
            IntVec<D> stdImage = minImage;
            mesh.shift(stdImage);
            std::cout << std::endl
                      << stdImage
                      << minImage
                      << Dbl(kSq[iw]) << Dbl(weights_[iw]*sum);
            #endif
         }
         UTIL_CHECK(std::abs(tot - 1.0) < 1.0E-8);
      }

   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D>::sample(long iStep)
   {
      // Preconditions
      UTIL_CHECK(isInitialized_);
      UTIL_CHECK(system().w().hasData());

      if (AnalyzerT::isAtInterval(iStep)) {

         // Compute W_{-}(r)
         RField<D> const & wa = system().w().rgrid(0);
         RField<D> const & wb = system().w().rgrid(1);
         //VecOp::subVV(wm_, wa, wb);
         //VecOp::mulEqS(wm_, 0.5);
         for (int i = 0; i < kSize_; ++i) {
             wm_[i] = 0.5 * ( wa[i] - wb[i] );
         }

         // Fourier transform W_{-}(r)
         system().domain().fft().forwardTransform(wm_, wk_);

         // Initialize bunch average values to zero
         int ib;
         for (ib = 0; ib < nBunch_; ++ib) {
            values_[ib] = 0.0;
         }

         const double vSystem  = system().domain().unitCell().volume();
         const double vMonomer = system().mixture().vMonomer();
         double chi = system().interaction().chi(0,1);
         double a_ = vSystem / (chi * chi * vMonomer * vMonomer);
         double b_ = 0.5 / (chi * vMonomer);

         // Compute structure factors, add to values_ array
         double value;
         for (int iw = 0; iw < kSize_; iw++) {
            value = a_ * absSq( wk_[iw] );
            value -= b_;
            ib = bunchIds_[iw];
            values_[ib] += weights_[iw] * value;
         }

         // Pass spherically averaged values to accumulators
         for (ib = 0; ib < nBunch_; ++ib) {
            accumulators_[ib].sample(values_[ib]);
         }

      }
   }

   /*
   * Output final results to output file.
   */
   template <int D>
   void BinaryStructureFactor<D>::output()
   {
      std::string const & outputFileName = AnalyzerT::outputFileName();
      AnalyzerT::fileMaster().openOutputFile(outputFileName, outputFile_);
      outputFile_ << "\t" << "q";
      outputFile_ << "\t" <<"\t" <<"S(q)";
      for (int i = 0; i < nBunch_; ++i) {
         outputFile_ << std::endl;
         outputFile_ << Dbl(wavenumbers_[i], 18, 8);
         outputFile_ << Dbl(accumulators_[i].average(), 18, 8);
      }
      outputFile_.close();
   }

}
}
#endif
