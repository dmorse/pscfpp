#ifndef RPC_FORCE_BIAS_MOVE_TPP
#define RPC_FORCE_BIAS_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMove.h"
#include "McMove.h" 
#include <rpc/simulate/mcmove/McSimulator.h>
#include <rpc/compressor/Compressor.h>
#include <rpc/System.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>


namespace Pscf {
namespace Rpc {

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
      system().compute();
      mcSimulator().computeWc();
      mcSimulator().computeCc();
      mcSimulator().computeDc();

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
      UTIL_CHECK(system().hasCFields());
      UTIL_CHECK(mcSimulator().hasWc());
      UTIL_CHECK(mcSimulator().hasCc());
      UTIL_CHECK(mcSimulator().hasDc());
      UTIL_CHECK(mcSimulator().hasHamiltonian());

      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Get current Hamiltonian
      double oldHamiltonian = mcSimulator().hamiltonian();

      // Save current state
      mcSimulator().saveMcState();

      // Clear both eigen-components of the fields and hamiltonian
      mcSimulator().clearData();

      attemptMoveTimer_.start();
      
      // Copy current W fields from parent system into wc_
      for (i = 0; i < nMonomer; ++i) {
         w_[i] = system().w().rgrid(i);
      }

      // Copy current derivative fields from into member variable dc_
      for (i = 0; i < nMonomer - 1; ++i) {
         dc_[i] = mcSimulator().dc(i);
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double vNode = vSystem/double(meshSize);
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_/vNode);

      // Modify local variables dwc_ and wc_
      // Loop over eigenvectors of projected chi matrix
      double dwd, dwr, evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = dc_[j];
         RField<D> & dwc = dwc_[j];
         for (k = 0; k < meshSize; ++k) {
            dwd = a*dc[k];
            dwr = b*random().gaussian();
            dwc[k] = dwd + dwr;
         }
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> const & dwc = dwc_[j];
            RField<D> & wn = w_[i];
            evec = mcSimulator().chiEvecs(j,i);
            for (k = 0; k < meshSize; ++k) {
               wn[k] += evec*dwc[k];
            }
         }
      }

      // Set modified fields in parent system
      system().setWRGrid(w_);
      mcSimulator().clearData();

      attemptMoveTimer_.stop();

      // Call compressor
      compressorTimer_.start();
      int compress = system().compressor().compress();
      UTIL_CHECK(compress == 0);
      compressorTimer_.stop();

      // Compute eigenvector components of current fields
      computeWcTimer_.start();
      mcSimulator().computeWc();
      mcSimulator().computeCc();
      mcSimulator().computeDc();
      computeWcTimer_.stop();

      // Evaluate new Hamiltonian
      computeHamiltonianTimer_.start();
      mcSimulator().computeHamiltonian();
      double newHamiltonian = mcSimulator().hamiltonian();
      double dH = newHamiltonian - oldHamiltonian;

      // Compute force bias 
      double dp, dm;
      double bias = 0.0;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & di = dc_[j];
         RField<D> const & df = mcSimulator().dc(j);
         RField<D> const & dwc = dwc_[j];
         for (k=0; k < meshSize; ++k) {
            dp = 0.5*(di[k] + df[k]);
            dm = 0.5*(di[k] - df[k]);
            bias += dp*( dwc[k] + mobility_*dm );
         }
      }
      bias *= vNode;
      computeHamiltonianTimer_.stop();

      // Accept or reject move
      decisionTimer_.start();
      bool accept = false;
      double weight = exp(bias - dH);
      accept = random().metropolis(weight);
      if (accept) {
          incrementNAccept();
          mcSimulator().clearMcState();
      } else {
          mcSimulator().restoreMcState();
          system().compute();
          mcSimulator().computeCc();
          mcSimulator().computeDc();
      }
      decisionTimer_.stop();
      totalTimer_.stop();

      return accept;
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
      out << "Real Move times contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif
