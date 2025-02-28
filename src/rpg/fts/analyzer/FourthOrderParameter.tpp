#ifndef RPG_FOURTH_ORDER_PARAMETER_TPP
#define RPG_FOURTH_ORDER_PARAMETER_TPP

#include "FourthOrderParameter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/System.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/cuda/resources.h>

#include <pscf/inter/Interaction.h>
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
#include <numeric>
#include <cmath>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D>::FourthOrderParameter(Simulator<D>& simulator, System<D>& system)
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      kSize_(1),
      hasAverage_(true),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("FourthOrderParameter"); }


   /*
   * Read parameters from file, and allocate memory.
   */
   template <int D>
   void FourthOrderParameter<D>::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      readOptional(in, "hasAverage", hasAverage_);
      readOptional(in,"nSamplePerBlock", nSamplePerBlock_);

      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      outputFile_ << "    chi       " << "FourthOrderParameter" << "\n";
   }

   /*
   * FourthOrderParameter setup
   */
   template <int D>
   void FourthOrderParameter<D>::setup()
   {
      //Check if the system is AB diblock copolymer
      const int nMonomer = system().mixture().nMonomer();
      if (nMonomer != 2) {
         UTIL_THROW("The FourthOrderParameter Analyzer is designed specifically for diblock copolymer system. Please verify the number of monomer types in your system.");
      }

      IntVec<D> const & dimensions = system().mesh().dimensions();

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

      // Allocate variables
      if (!isInitialized_){
         wc0_.allocate(dimensions);
         wK_.allocate(dimensions);
         prefactor_.allocate(kMeshDimensions_);
         VecOp::eqS(prefactor_, 0);
      }

      isInitialized_ = true;

      // Clear accumulators
      if (hasAverage_){
         accumulator_.clear();
      }

      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }

      computePrefactor();
   }

   /*
   * Increment structure factors for all wavevectors and modes.
   */
   template <int D>
   void FourthOrderParameter<D>::sample(long iStep)
   {
      if (!isAtInterval(iStep)) return;
      computeFourthOrderParameter();

      if (hasAverage_){
         accumulator_.sample(FourthOrderParameter_);
      }

      double chi =  system().interaction().chi(0,1);
      UTIL_CHECK(outputFile_.is_open());
      outputFile_ << Dbl(chi);
      outputFile_ << Dbl(FourthOrderParameter_);
      outputFile_<< "\n";
   }

   template <int D>
   void FourthOrderParameter<D>::computeFourthOrderParameter()
   {
      UTIL_CHECK(system().w().hasData());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      const int meshSize = system().domain().mesh().size();
      RField<D> psi;
      psi.allocate(kMeshDimensions_);

      // Convert W_(r) to fourier mode W_(k)
      VecOp::eqV(wc0_, simulator().wc(0));
      system().fft().forwardTransform(wc0_, wK_);

      // psi = |wK_|^4
      VecOp::sqSqNormV(psi, wK_);

      // W_(k)^4 * weight factor
      VecOp::mulEqV(psi, prefactor_);

      // Get sum over all wavevectors
      FourthOrderParameter_ = Reduce::sum(psi);
      FourthOrderParameter_ = std::pow(FourthOrderParameter_, 0.25);
   }

   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
   {
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      HostDArray<cudaReal> prefactor_h(kSize_);
      for (int i = 0; i < kSize_; ++i){
         prefactor_h[i] = 0;
      }
      IntVec<D> G;
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      DArray<IntVec<D>> GminList;
      GminList.allocate(kSize_);
      MeshIterator<D> itr(kMeshDimensions_);
      MeshIterator<D> searchItr(kMeshDimensions_);

      // Calculate GminList
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         GminList[itr.rank()] = Gmin;
      }

      // Compute weight factor for each G wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If the weight factor of the current wavevector has not been assigned
         if (prefactor_h[itr.rank()] == 0){
            Gmin = GminList[itr.rank()];

            // Compute inverse of wavevector
            nGmin.negate(Gmin);

            // Search for inverse of wavevector
            searchItr = itr;
            for (; !searchItr.atEnd(); ++searchItr){
               if (nGmin == GminList[searchItr.rank()]){
                  prefactor_h[itr.rank()] = 1.0/2.0;
                  prefactor_h[searchItr.rank()] = 1.0/2.0;
                  inverseFound = true;
               }
            }

            if (inverseFound == false){
               prefactor_h[itr.rank()]  = 1.0;
            }

         }

      }

      // Copy the weight factor from cpu(host) to gpu(device)
      prefactor_ = prefactor_h;
   }

   /*
   * Output final results to output file.
   */
   template <int D>
   void FourthOrderParameter<D>::output()
   {
      if (hasAverage_){
         Log::file() << std::endl;
         Log::file() << "At chi = " << system().interaction().chi(0,1) << "\n";
         Log::file() << "Time average of the FourthOrderParameter is: "
                     << Dbl(accumulator_.average())
                     << " +- " << Dbl(accumulator_.blockingError(), 9, 2)
                     << "\n";
      }

   }

}
}
#endif
