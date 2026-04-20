#ifndef RPC_STRUCTURE_FACTOR_BUNCHED_TPP
#define RPC_STRUCTURE_FACTOR_BUNCHED_TPP

#include "StructureFactorBunched.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/FFT.h>
#include <prdc/cpu/WaveList.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/misc/FileMaster.h>
//#include <util/misc/ioUtil.h>
#include <util/format/Dbl.h>
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
   StructureFactorBunched<D>::StructureFactorBunched(Simulator<D>& simulator,
                                         System<D>& system)
    : Analyzer<D>(simulator, system),
      kMeshDimensions_(0),
      kSize_(0),
      nBunch_(0),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  ParamComposite::setClassName("StructureFactorBunched"); }

   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void StructureFactorBunched<D>::readParameters(std::istream& in)
   {
      // Precondition: Require that the system has two monomer types
      UTIL_CHECK(system().mixture().nMonomer() == 2);

      AnalyzerT::readInterval(in);
      AnalyzerT::readOutputFileName(in);
      ParamComposite::readOptional(in,"nSamplePerBlock", nSamplePerBlock_);
      isInitialized_ = true;
   }

   /*
   * StructureFactorBunched setup
   */
   template <int D>
   void StructureFactorBunched<D>::setup()
   {
      UTIL_CHECK(isInitialized_);

      // Compute Fourier space kMeshDimensions_ and kSize_
      IntVec<D> const & rMeshDimensions = system().domain().mesh().dimensions();
      FFT<D>::computeKMesh(rMeshDimensions, kMeshDimensions_, kSize_);

      // If needed, allocate arrays indexed by wave id
      if (!wm_.isAllocated()){
         wm_.allocate(rMeshDimensions);
         wk_.allocate(rMeshDimensions);
         bunchIds_.allocate(kSize_);
         weights_.allocate(kSize_);
         values_.allocate(kSize_);
      }
      UTIL_CHECK(wm_.capacity() == kSize_);
      UTIL_CHECK(wk_.capacity() == kSize_);
      UTIL_CHECK(bunchIds_.capacity() == kSize_);
      UTIL_CHECK(weights_.capacity() == kSize_);
      UTIL_CHECK(values_.capacity() == kSize_);

      // Sort waves and set nBunch
      system().waveList().sortWaves();
      nBunch_ = system().waveList().nBunch();

      // If needed, allocate arrays indexed by bunch id
      if (accumulators_.isAllocated()) {
         UTIL_CHECK(accumulators_.capacity() == wavenumbers_.capacity());
         if (accumulators_.capacity() < nBunch_) {
            accumulators_.deallocate();
            wavenumbers_.deallocate();
         }
      }
      if (!accumulators_.isAllocated()) {
         UTIL_CHECK(!wavenumbers_.isAllocated());
         accumulators_.allocate(nBunch_);
         wavenumbers_.allocate(nBunch_);
      }
      int n = accumulators_.capacity();
      UTIL_CHECK(wavenumbers_.capacity() == n);
      UTIL_CHECK(n >= nBunch_);

      // Initialize accumulators
      int ib;
      for (ib = 0; ib < n; ++ib) {
         accumulators_[ib].setNSamplePerBlock(nSamplePerBlock_);
         accumulators_[ib].clear();
         wavenumbers_[ib] = 0.0;
      }

      // Initialize bunchIds_, weights_, and wavenumbers_
      WaveList<D> const & waveList = system().domain().waveList();
      Array< double > const & kSq = waveList.kSq();
      Array< int > const & sortedIds = waveList.sortedIds();
      GArray< Pair<int> > const & bunches = waveList.sortedBunches();
      Array<bool> const & implicit = waveList.implicitInverse();
      double wavenumberSq, count, sum, tot;
      int begin, end, iw, k;
      for (ib = 0; ib < nBunch_; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];
         wavenumberSq = kSq[begin];
         UTIL_CHECK(wavenumberSq >= 0.0);       
         wavenumbers_[ib] = std::sqrt(wavenumberSq);
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
         tot = 0.0;
         for (k = begin; k < end; ++k) {
            iw = sortedIds[k];
            weights_[iw] /= sum;
            tot += weights_[iw];
         }
         UTIL_CHECK(std::abs(tot - 1.0) < 1.0E-6);
      }

   }

   /*
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void StructureFactorBunched<D>::sample(long iStep)
   {
      // Preconditions
      UTIL_CHECK(isInitialized_);
      UTIL_CHECK(system().w().hasData());

      if (AnalyzerT::isAtInterval(iStep)) {

         // Compute W_{-}(r)
         double wa, wb;
         for (int i = 0; i < kSize_; ++i) {
             wa = system().w().rgrid(0)[i];
             wb = system().w().rgrid(1)[i];
             wm_[i] = 0.5 * ( wa - wb );
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

         // Pass square magnitudes of Fourier components to accumulators
         double wr, wi, value;
         for (int iw = 0; iw < kSize_; iw++) {
            wr = wk_[iw][0];
            wi = wk_[iw][1];
            value = a_ * ( wr*wr + wi*wi ) - b_;
            ib = bunchIds_[iw];
            values_[ib] += weights_[iw] * value;
         }

         // Pass bunch average values to accumulators
         for (ib = 0; ib < nBunch_; ++ib) {
            accumulators_[ib].sample(values_[ib]);
         }

      }
   }

   #if 0
   template <int D>
   void StructureFactorBunched<D>::computeStructureFactor()
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
   #endif

   /*
   * Output final results to output file.
   */
   template <int D>
   void StructureFactorBunched<D>::output()
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
