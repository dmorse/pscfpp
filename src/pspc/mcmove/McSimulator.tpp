#ifndef PSPC_MC_SIMULATOR_TPP
#define PSPC_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"
#include <pspc/System.h>
#include <pspc/mcmove/McMoveFactory.h>

#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McSimulator<D>::McSimulator(System<D>& system)
   : Manager< McMove<D> >(),
     systemPtr_(&system),
     randomPtr_(&system.random())
   {  setClassName("McSimulator"); }

   /*
   * Destructor.
   */
   template <int D>
   McSimulator<D>::~McSimulator()
   {}

   /*
   * Return a pointer to a new McMoveFactory object.
   */
   template <int D>
   Factory< McMove<D> >* McSimulator<D>::newDefaultFactory() const
   {  return new McMoveFactory<D>(*systemPtr_); }

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void McSimulator<D>::readParameters(std::istream &in)
   {
      Manager< McMove<D> >::readParameters(in);

      // Allocate and store probabilities
      probabilities_.allocate(size());
      double  totalProbability = 0.0;
      int     iMove;
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = (*this)[iMove].probability();
         totalProbability += probabilities_[iMove];
      }

      // Allocate and store and normalize probabilities
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = probabilities_[iMove]/totalProbability;
         (*this)[iMove].setProbability(probabilities_[iMove]);
      }
   }

   /*
   * Initialize all moves just prior to a run.
   */
   template <int D>
   void McSimulator<D>::setup()
   {
      for (int iMove = 0; iMove < size(); ++iMove) {
         (*this)[iMove].setup();
      }
   }

   /*
   * Choose a McMove at random.
   */
   template <int D>
   McMove<D>& McSimulator<D>::chooseMove()
   {
      int iMove;
      iMove = randomPtr_->drawFrom(&probabilities_[0], size());
      return (*this)[iMove];
   }

   /*
   * Output statistics for every move.
   */
   template <int D>
   void McSimulator<D>::output()
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].output();
      }
   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void McSimulator<D>::simulate(int nStep)
   {

      setup();
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for (int iStep = 0; iStep < nStep; ++iStep) {

         // Choose and attempt an McMove
         chooseMove().move();

      }
      timer.stop();
      double time = timer.time();

      // Output results of move statistics to files
      output();

      // Output time for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep         " << nStep << std::endl;
      Log::file() << "run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep  " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << std::endl;

      // Print McMove acceptance statistics
      long attempt;
      long accept;
      using namespace std;
      Log::file() << "Move Statistics:" << endl << endl;
      Log::file() << setw(32) << left <<  "Move Name"
           << setw(12) << right << "Attempted"
           << setw(12) << right << "Accepted"
           << setw(15) << right << "AcceptRate"
           << endl;
      int nMove = size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = (*this)[iMove].nAttempt();
         accept  = (*this)[iMove].nAccept();
         Log::file() << setw(32) << left
              << (*this)[iMove].className()
              << setw(12) << right << attempt
              << setw(12) << accept
              << setw(15) << fixed << setprecision(6)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

   }

   /*
   * Save the current Monte-Carlo state.
   *
   * Used before attempting a Monte-Carlo move.
   */
   template <int D>
   void McSimulator<D>::saveMcState()
   {
      //UTIL_CHECK(system().isAllocatedRGrid());
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(!mcState_.hasData);

      int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         mcState_.w[i] = system().w().rgrid(i);
      }
      mcState_.mcHamiltonian  = mcHamiltonian_;
      mcState_.hasData = true;
   }

   /*
   * Restore a saved Monte-Carlo state.
   *
   * Used when an attempted Monte-Carlo move is rejected.
   */
   template <int D>
   void McSimulator<D>::restoreMcState()
   {
      system().setWRGrid(mcState_.w); 
      mcHamiltonian_ = mcState_.mcHamiltonian;
      mcState_.hasData = false;
      hasMcHamiltonian_ = true;
   }

   /*
   * Compute Monte Carlo Hamiltonian.
   */
   template <int D>
   void McSimulator<D>::computeMcHamiltonian()
   {
      UTIL_CHECK(system().domain().basis().isInitialized());
      UTIL_CHECK(system().w().hasData());
      //UTIL_CHECK(hasCFields_);

      system().computeFreeEnergy();

      /// Computation, sets mcHamiltonian_

      hasMcHamiltonian_ = true;
   }

}
}
#endif
