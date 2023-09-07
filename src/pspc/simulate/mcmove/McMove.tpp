#ifndef PSPC_MC_MOVE_TPP
#define PSPC_MC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"

#include <pspc/System.h>
#include <util/archives/Serializable_includes.h>
#include <pspc/compressor/Compressor.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMove<D>::McMove(McSimulator<D>& mcSimulator)
    : mcSimulatorPtr_(&mcSimulator),
      systemPtr_(&(mcSimulator.system())),
      randomPtr_(&(mcSimulator.random())),
      nAttempt_(0),
      nAccept_(0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   McMove<D>::~McMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void McMove<D>::readParameters(std::istream &in)
   {}

   /*
   * Read the probability from file.
   */
   template <int D>
   void McMove<D>::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }

   /*
   * Setup at beginning of loop.
   *
   * Trivial default implementation - initializes counters.
   */
   template <int D>
   void McMove<D>::setup()
   {
      nAttempt_ = 0;
      nAccept_  = 0;
      clearTimers();
   }

   /*
   * Trivial default implementation - always returns false.
   */
   template <int D>
   bool McMove<D>::move()
   {
      totalTimer_.start();
      incrementNAttempt();

      // Get current Hamiltonian
      double oldHamiltonian = mcSimulator().mcHamiltonian();

      // Save current state
      mcSimulator().saveMcState();

      // Clear both eigen-components of the fields and mcHamiltonian
      mcSimulator().clearData();
      
      // Attempt modification
      attemptMoveTimer_.start();
      attemptMove();
      attemptMoveTimer_.stop();
      
      // Call compressor
      compressorTimer_.start();
      system().compressor().compress();
      compressorTimer_.stop();
      
      // Compute eigenvector components of the current w fields
      computeWCTimer_.start();
      mcSimulator().computeWC();
      computeWCTimer_.stop();
      
      // Evaluate new Hamiltonian
      computeMcHamiltonianTimer_.start();
      mcSimulator().computeMcHamiltonian();
      double newHamiltonian = mcSimulator().mcHamiltonian();
      computeMcHamiltonianTimer_.stop();

      // Accept or reject move
      decisionTimer_.start();
      bool accept = false;
      double weight = exp(-(newHamiltonian - oldHamiltonian));
      accept = random().metropolis(weight);
      if (accept) {
          incrementNAccept();
          mcSimulator().clearMcState();
      } else {
          mcSimulator().restoreMcState();
      }
      decisionTimer_.stop();
      totalTimer_.stop();
      
      return accept;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void McMove<D>::output()
   {}
   
   template<int D>
   void McMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "McMove times contributions:\n";
      // Output timing results, if requested.
      double total = totalTimer_.time();
      out << "Attempt Move:         "
          << Dbl(attemptMoveTimer_.time(), 9, 3)  << " s,  "
          << Dbl(attemptMoveTimer_.time()/total, 9, 3) << " %" << "\n";
      out << "Compressor: "
          << Dbl(compressorTimer_.time(), 9, 3)  << " s,  "
          << Dbl(compressorTimer_.time()/total, 9, 3) << " %" << "\n";
      out << "Compute eigen-components of the fields:  "
          << Dbl(computeWCTimer_.time(), 9, 3)  << " s,  "
          << Dbl(computeWCTimer_.time()/total, 9, 3) << " %" << "\n";
      out << "Compute Hamiltonian: "
          << Dbl(computeMcHamiltonianTimer_.time(), 9, 3)  << " s,  "
          << Dbl(computeMcHamiltonianTimer_.time()/total, 9, 3) << " %" << "\n";
      out << "Make decision (accept or reject):       "
          << Dbl(decisionTimer_.time(), 9, 3)  << " s,  "
          << Dbl(decisionTimer_.time()/total, 9, 3) << " %" << "\n";
      out << "total time:           "
          << Dbl(total, 9, 3) << " s  \n";
      out << "\n";
   }
   
   template<int D>
   void McMove<D>::clearTimers()
   {
      attemptMoveTimer_.clear();
      compressorTimer_.clear();
      computeWCTimer_.clear();
      computeMcHamiltonianTimer_.clear();
      decisionTimer_.clear();
      totalTimer_.clear();
   }
      

}
}
#endif
