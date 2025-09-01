#ifndef RPG_MC_MOVE_TPP
#define RPG_MC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"

#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>
#include <pscf/cuda/CudaRandom.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMove<D>::McMove(McSimulator<D>& simulator)
    : simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      randomPtr_(&(simulator.random())),
      cudaRandomPtr_(&(simulator.cudaRandom())),
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
      nFail_ = 0;
      clearTimers();
      simulator().computeWc();
      if (simulator().needsCc() || simulator().needsDc()){
         system().compute();
         simulator().computeCc();
      }
      if (simulator().needsDc()){
         simulator().computeDc();
      }
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
      double oldHamiltonian = simulator().hamiltonian();

      // Save current state
      simulator().saveState();

      // Clear both eigen-components of the fields and hamiltonian
      simulator().clearData();
      
      // Attempt modification
      attemptMoveTimer_.start();
      attemptMove();
      attemptMoveTimer_.stop();
      
      // Call compressor
      compressorTimer_.start();
      int compress = simulator().compressor().compress();
      UTIL_CHECK(system().c().hasData());
      compressorTimer_.stop();
      
      bool isConverged = false;
      if (compress != 0){
         incrementNFail();
         simulator().restoreState();
      } else {
         isConverged = true;
         
         // Compute eigenvector components of the current w fields
         componentTimer_.start();
         simulator().computeWc();
         // Compute cc fields if any move require cc fields
         if (simulator().needsCc() || simulator().needsDc()){
            //system().compute();
            UTIL_CHECK(system().c().hasData());
            simulator().computeCc();
         }
         // Compute dc fields if any move require dc fields
         if (simulator().needsDc()){
            simulator().computeDc();
         }
         componentTimer_.stop();
      
         // Evaluate new Hamiltonian
         hamiltonianTimer_.start();
         simulator().computeHamiltonian();
         double newHamiltonian = simulator().hamiltonian();
         hamiltonianTimer_.stop();
         
         // Accept or reject move
         bool accept = false;
         decisionTimer_.start();
         double weight = exp(-(newHamiltonian - oldHamiltonian));
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
   void McMove<D>::output()
   {}
   
   template<int D>
   void McMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      double total = totalTimer_.time();
      out << "                          "
          << "Total" << std::setw(17) << "Per Move" << std::setw(14) << "Fraction" << "\n";
      out << "Attempt Move:             "
          << Dbl(attemptMoveTimer_.time(), 9, 3)  << " s,  "
          << Dbl(attemptMoveTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(attemptMoveTimer_.time()/total, 9, 3) << "\n";
      out << "Compressor:               "
          << Dbl(compressorTimer_.time(), 9, 3)  << " s,  "
          << Dbl(compressorTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(compressorTimer_.time()/total, 9, 3) << "\n";
      out << "Compute eigen-components: "
          << Dbl(componentTimer_.time(), 9, 3)  << " s,  "
          << Dbl(componentTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(componentTimer_.time()/total, 9, 3) << "\n";
      out << "Compute Hamiltonian:      "
          << Dbl(hamiltonianTimer_.time(), 9, 3)  << " s,  "
          << Dbl(hamiltonianTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(hamiltonianTimer_.time()/total, 9, 3) << "\n";
      out << "Accept or Reject:         "
          << Dbl(decisionTimer_.time(), 9, 3)  << " s,  "
          << Dbl(decisionTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(decisionTimer_.time()/total, 9, 3) << "\n";
      out << "total time:               "
          << Dbl(total, 9, 3) << " s,  "
          << Dbl(total/nAttempt_, 9, 3) << " s  \n";
      out << "\n";
   }
   
   template<int D>
   void McMove<D>::clearTimers()
   {
      attemptMoveTimer_.clear();
      compressorTimer_.clear();
      componentTimer_.clear();
      hamiltonianTimer_.clear();
      decisionTimer_.clear();
      totalTimer_.clear();
   }
      

}
}
#endif
