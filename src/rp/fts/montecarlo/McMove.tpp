#ifndef RP_MC_MOVE_TPP
#define RP_MC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"

#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/simulator/SimState.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/system/System.h>
#include <rp/field/CFields.h>
#include <prdc/field/RField.h>
#include <util/random/Random.h>
#include <util/archives/Serializable_includes.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   McMove<D,T>::McMove(McSimulator<D,T>& simulator)
    : simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      randomPtr_(&(simulator.random())),
      vecRandomPtr_(&(simulator.vecRandom())),
      probability_(0.0),
      nAttempt_(0),
      nAccept_(0),
      nFail_(0)
   {}

   /*
   * Read parameters - empty default implementation.
   */
   template <int D, class T>
   void McMove<D,T>::readParameters(std::istream &in)
   {}

   /*
   * Read the probability from file.
   */
   template <int D, class T>
   void McMove<D,T>::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }

   /*
   * Setup at beginning of loop.
   *
   * Trivial default implementation - initializes counters.
   */
   template <int D, class T>
   void McMove<D,T>::setup()
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
   * Standard Metropolis MC move.
   */
   template <int D, class T>
   bool McMove<D,T>::move()
   {
      // Start timers
      totalTimer_.start();
      attemptMoveTimer_.start();
      incrementNAttempt();

      // Get current Hamiltonian
      double oldHamiltonian = simulator().hamiltonian();
   
      // Save current state
      simulator().saveState();

      // Clear both eigen-components of the fields and hamiltonian
      simulator().clearData();

      // Attempt modification of exchange fields
      attemptMove();
      attemptMoveTimer_.stop();

      // Call compressor to modify pressure-like fields
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

         // Compute eigenvector components of the current fields
         componentTimer_.start();
         simulator().computeWc();
         // Compute cc fields if any move require cc fields
         if (simulator().needsCc() || simulator().needsDc()){
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
         decisionTimer_.start();
         bool accept = false;
         double weight = exp(-(newHamiltonian - oldHamiltonian));
         if (newHamiltonian - oldHamiltonian <= 1e-10){
            accept = true;
         } else{
            accept = random().metropolis(weight);
         }
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
   template <int D, class T>
   void McMove<D,T>::output()
   {}

   /*
   * Output timing information at tend of simulation.
   */
   template<int D, class T>
   void McMove<D,T>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << className();
      out << " time contributions:\n";
      outputTimerData(out);
   }

   /*
   * Output raw timing information, with class name label.
   */
   template<int D, class T>
   void McMove<D,T>::outputTimerData(std::ostream& out)
   {
      // Output timing results, if requested.
      double total = totalTimer_.time();
      out << "                          " << "Total" 
          << std::setw(17) << "Per Move" 
          << std::setw(14) << "Fraction" << "\n";
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

   /*
   * Clear all timers.
   */
   template<int D, class T>
   void McMove<D,T>::clearTimers()
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
