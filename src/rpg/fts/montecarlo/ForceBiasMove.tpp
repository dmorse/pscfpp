#ifndef RPG_FORCE_BIAS_MOVE_TPP
#define RPG_FORCE_BIAS_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMove.h"
#include "McMove.h" 
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/fts/VecOpFts.h>
#include <rpg/System.h>
#include <prdc/cuda/resources.h>
#include <pscf/cuda/CudaRandom.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D>::ForceBiasMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      w_(),
      dwc_(),
      mobility_(0.0)
   {  setClassName("ForceBiasMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   ForceBiasMove<D>::~ForceBiasMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void ForceBiasMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);

      // Read Brownian dynamics mobility parameter
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dc_.allocate(nMonomer-1);
      dwc_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         dc_[i].allocate(meshDimensions);
         dwc_[i].allocate(meshDimensions);
      }

   }
   
   template <int D>
   void ForceBiasMove<D>::setup()
   {  
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      IntVec<D> dimensions = system().domain().mesh().dimensions();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dc_.capacity() == nMonomer-1);
      UTIL_CHECK(dwc_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dc_[i].capacity() == meshSize);
         UTIL_CHECK(dwc_[i].capacity() == meshSize);
      }

      McMove<D>::setup();
      biasField_.allocate(dimensions);
      gaussianField_.allocate(dimensions);
   }
 
   /*
   * Attempt unconstrained move
   */
   template <int D>
   bool ForceBiasMove<D>::move()
   {
      totalTimer_.start();
      incrementNAttempt();

      // Preconditions
      UTIL_CHECK(simulator().hasWc());
      UTIL_CHECK(simulator().hasDc());
      UTIL_CHECK(simulator().hasHamiltonian());

      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      

      // Get current Hamiltonian
      double oldHamiltonian = simulator().hamiltonian();

      // Save current state
      simulator().saveState();

      // Clear both eigen-components of the fields and hamiltonian
      simulator().clearData();

      attemptMoveTimer_.start();
      
      // Copy current W fields from parent system into wc_
      for (i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Copy current derivative fields from into member variable dc_
      for (i = 0; i < nMonomer - 1; ++i) {
         VecOp::eqV(dc_[i], simulator().dc(i));
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double vNode = vSystem/double(meshSize);
      double a = -1.0*mobility_;
      double b = sqrt(2.0*mobility_/vNode);
      
      // Constants for normal distribution
      double stddev = 1.0;
      double mean = 0;

      // Modify local variables dwc_ and wc_
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> & dwc = dwc_[j];
         
         // Generate normal distributed random floating point numbers
         cudaRandom().normal(gaussianField_, stddev, mean);
         
         // dwc
         VecOp::addVcVc(dwc, dc_[j], a, gaussianField_, b);
         
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> const & dwc = dwc_[j];
            RField<D> & wn = w_[i];
            evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(wn, dwc, evec);
         }
         
      }

      // Set modified fields in parent system
      system().setWRGrid(w_);
      simulator().clearData();

      attemptMoveTimer_.stop();

      // Call compressor
      compressorTimer_.start();
      int compress = simulator().compressor().compress();
      UTIL_CHECK(system().hasCFields());
      compressorTimer_.stop();
      
      bool isConverged = false;
      if (compress != 0){
         incrementNFail();
         simulator().restoreState();
      } else {
         isConverged = true;
         
         // Compute eigenvector components of current fields
         componentTimer_.start();
         simulator().computeWc();
         UTIL_CHECK(system().hasCFields());
         simulator().computeCc();
         simulator().computeDc();
         componentTimer_.stop();

         // Evaluate new Hamiltonian
         hamiltonianTimer_.start();
         simulator().computeHamiltonian();
         double newHamiltonian = simulator().hamiltonian();
         double dH = newHamiltonian - oldHamiltonian;

         // Compute force bias 
         double bias = 0.0;
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D> const & di = dc_[j];
            RField<D> const & df = simulator().dc(j);
            RField<D> const & dwc = dwc_[j];
         
            VecOpFts::computeForceBias(biasField_, di, df, dwc, mobility_);
            bias += Reduce::sum(biasField_);
         }
         bias *= vNode;
         hamiltonianTimer_.stop();

         // Accept or reject move
         decisionTimer_.start();
         bool accept = false;
         double weight = exp(bias - dH);
         accept = random().metropolis(weight);
         if (accept) {
            incrementNAccept();
            simulator().clearState();
         } else {
            simulator().restoreState();
         }
         decisionTimer_.stop();
      }
      
      totalTimer_.stop();
      return isConverged;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void ForceBiasMove<D>::output()
   {}

   template<int D>
   void ForceBiasMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "ForceBiasMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif
