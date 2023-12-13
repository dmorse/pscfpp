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
#include <pspc/simulate/mcmove/McMoveFactory.h>
#include <pspc/simulate/analyzer/AnalyzerFactory.h>
#include <pspc/simulate/trajectory/TrajectoryReader.h>
#include <pspc/simulate/trajectory/TrajectoryReaderFactory.h>
#include <pspc/compressor/Compressor.h>

#include <util/random/Random.h>
#include <util/misc/Timer.h>
#include <util/global.h>

#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McSimulator<D>::McSimulator(System<D>& system)
    : Simulator<D>(system),
      mcMoveManager_(*this, system),
      analyzerManager_(*this, system),
      trajectoryReaderFactoryPtr_(0)
   { 
      setClassName("McSimulator");   
      trajectoryReaderFactoryPtr_ 
             = new TrajectoryReaderFactory<D>(system);
   }

   /*
   * Destructor.
   */
   template <int D>
   McSimulator<D>::~McSimulator()
   {
      if (trajectoryReaderFactoryPtr_) {
         delete trajectoryReaderFactoryPtr_;
      }
   }

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void McSimulator<D>::readParameters(std::istream &in)
   {
      // Initialize Simulator<D> base class
      allocate();
      analyzeChi();

      // Read block of mc move data inside
      readParamComposite(in, mcMoveManager_);

      // Read block of analyzer
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);

   }

   /*
   * Initialize just prior to a run.
   */
   template <int D>
   void McSimulator<D>::setup()
   {  
      UTIL_CHECK(system().w().hasData());

      // Allocate mcState_, if necessary.
      if (!mcState_.isAllocated) {
         const int nMonomer = system().mixture().nMonomer();
         const IntVec<D> dimensions = system().domain().mesh().dimensions();
         mcState_.allocate(nMonomer, dimensions);
      }
   
      // Eigenanalysis of the projected chi matrix.
      analyzeChi();

      // Compute field components and MC Hamiltonian for initial state
      system().compute();
      computeWc();
      computeHamiltonian();
      mcMoveManager_.setup();
      if (analyzerManager_.size() > 0){
         analyzerManager_.setup();
      }
      
   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void McSimulator<D>::simulate(int nStep)
   {
      UTIL_CHECK(mcMoveManager_.size() > 0);

      setup();
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      Timer analyzerTimer;
      timer.start();
      for (iStep_ = 0; iStep_ < nStep; ++iStep_) {

         // Analysis (if any)
         analyzerTimer.start();
         if (Analyzer<D>::baseInterval != 0) {
            if (iStep_ % Analyzer<D>::baseInterval == 0) {
               if (analyzerManager_.size() > 0) {
                  analyzerManager_.sample(iStep_);
               }
            }
         }
         analyzerTimer.stop();

         // Choose and attempt an McMove
         mcMoveManager_.chooseMove().move();

      }

      // Analysis (if any)
      analyzerTimer.start();
      if (Analyzer<D>::baseInterval != 0) {
         if (iStep_ % Analyzer<D>::baseInterval == 0) {
            if (analyzerManager_.size() > 0) {
               analyzerManager_.sample(iStep_);
            }
         }
      }
      analyzerTimer.stop();

      timer.stop();
      double time = timer.time();
      double analyzerTime = analyzerTimer.time();

      // Output results of move statistics to files
      mcMoveManager_.output();
      if (Analyzer<D>::baseInterval > 0){
         analyzerManager_.output();
      }

      // Output number of times MDE has been solved for the simulation run
      Log::file() << std::endl;
      Log::file() << "MDE counter   " 
                  << system().compressor().mdeCounter() << std::endl;
      Log::file() << std::endl;
      
      // Output times for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep               " << nStep << std::endl;
      Log::file() << "Total run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep        " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << "Analyzer run time   " << analyzerTime
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
      int nMove = mcMoveManager_.size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = mcMoveManager_[iMove].nAttempt();
         accept  = mcMoveManager_[iMove].nAccept();
         Log::file() << setw(32) << left
              << mcMoveManager_[iMove].className()
              << setw(12) << right << attempt
              << setw(12) << accept
              << setw(15) << fixed << setprecision(6)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

   }

   /*
   * Save the current Monte-Carlo state prior to an attempted MC move.
   *
   * Invoked before each attempted Monte-Carlo move.
   */
   template <int D>
   void McSimulator<D>::saveMcState()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(hasWc());
      UTIL_CHECK(hasHamiltonian());
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(!mcState_.hasData);

      // Set fields
      int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         mcState_.w[i] = system().w().rgrid(i);
         mcState_.wc[i] = wc(i);
      }

      // Set Hamiltonian
      mcState_.hamiltonian  = hamiltonian();
      mcState_.idealHamiltonian  = idealHamiltonian();
      mcState_.fieldHamiltonian  = fieldHamiltonian();

      mcState_.hasData = true;
   }

   /*
   * Restore a saved Monte-Carlo state.
   *
   * Invoked after an attempted Monte-Carlo move is rejected.
   */
   template <int D>
   void McSimulator<D>::restoreMcState()
   {
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(mcState_.hasData);
      const int nMonomer = system().mixture().nMonomer();

      // Restore fields
      system().setWRGrid(mcState_.w); 
      for (int i = 0; i < nMonomer; ++i) {
         wc_[i] = mcState_.wc[i];
      }
      hasWc_ = true;

      // Restore Hamiltonian and components
      system().setWRGrid(mcState_.w); 
      hamiltonian_ = mcState_.hamiltonian;
      idealHamiltonian_ = mcState_.idealHamiltonian;
      fieldHamiltonian_ = mcState_.fieldHamiltonian;
      hasHamiltonian_ = true;

      mcState_.hasData = false;
   }
 
   /*
   * Clear the saved Monte-Carlo state.
   *
   * Invoked when an attempted Monte-Carlo move is accepted.
   */
   template <int D>
   void McSimulator<D>::clearMcState()
   {  mcState_.hasData = false; }

   /*
   * Open, read and analyze a trajectory file
   */
   template <int D>
   void McSimulator<D>::analyze(int min, int max,
                                std::string classname,
                                std::string filename)
   {
      // Preconditions
      UTIL_CHECK(min >= 0);
      UTIL_CHECK(max >= min);
      UTIL_CHECK(Analyzer<D>::baseInterval > 0);
      UTIL_CHECK(analyzerManager_.size() > 0);
      
      // Construct TrajectoryReader
      TrajectoryReader<D>* trajectoryReaderPtr;
      trajectoryReaderPtr = trajectoryReaderFactory().factory(classname);
      if (!trajectoryReaderPtr) {
         std::string message;
         message = "Invalid TrajectoryReader class name " + classname;
         UTIL_THROW(message.c_str());
      }

      // Open trajectory file
      trajectoryReaderPtr->open(filename);
      trajectoryReaderPtr->readHeader();

      // Main loop over trajectory frames
      Timer timer;
      bool hasFrame = true;
      timer.start();
      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         hasFrame = trajectoryReaderPtr->readFrame();
         if (hasFrame) {
            clearData();

            // Initialize analyzers 
            if (iStep_ == min) {
               analyzerManager_.setup();
            }

            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager_.sample(iStep_);
            }
         }
      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;
      int nFrames = iStep_ - min;
      trajectoryReaderPtr->close();
      delete trajectoryReaderPtr;

      // Output results of all analyzers to output files
      analyzerManager_.output();

      // Output number of frames and times
      Log::file() << std::endl;
      Log::file() << "# of frames   " << nFrames << std::endl;
      Log::file() << "run time      " << timer.time() 
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames) 
                  << "  sec" << std::endl;
      Log::file() << std::endl;

   }
   
   /*
   * Output McMoveManager timer results.
   */ 
   template<int D>
   void McSimulator<D>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << "McSimulator times contributions:\n";
      mcMoveManager_.outputTimers(out);
   }
  
   /*
   * Clear all McMoveManager timers.
   */ 
   template<int D>
   void McSimulator<D>::clearTimers()
   {  mcMoveManager_.clearTimers(); }

}
}
#endif
