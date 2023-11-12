#ifndef PSPG_MC_SIMULATOR_TPP
#define PSPG_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <pspg/System.h>
#include <pspg/simulate/mcmove/McMoveFactory.h>
#include <pspg/simulate/analyzer/AnalyzerFactory.h>
#include <pspg/simulate/trajectory/TrajectoryReader.h>
#include <pspg/simulate/trajectory/TrajectoryReaderFactory.h>
#include <pspg/compressor/Compressor.h>

//#include <util/random/Random.h>
#include <util/misc/Timer.h>
#include <util/global.h>

#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Pscf::Prdc::Cuda;

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
      trajectoryReaderFactoryPtr_ = new TrajectoryReaderFactory<D>(system);
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
      // Read block of mc move data inside
      readParamComposite(in, mcMoveManager_);

      // Read block of analyzer
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);

      Simulator<D>::readParameters(in);
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
         int meshSize = system().domain().mesh().size();
         mcState_.allocate(nMonomer, meshSize);
      }
   
      // Eigenanalysis of the projected chi matrix.
      analyzeChi();

      // Compute field components and MC Hamiltonian for initial state
      system().compute();
      computeWC();
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
      timer.start();
      for (int iStep = 0; iStep < nStep; ++iStep) {

         // Analysis (if any)
         if (Analyzer<D>::baseInterval != 0) {
            if (iStep % Analyzer<D>::baseInterval == 0) {
               if (analyzerManager_.size() > 0) {
                  iStep_ = iStep;
                  analyzerManager_.sample(iStep);
               }
            }
         }

         // Choose and attempt an McMove
         mcMoveManager_.chooseMove().move();

      }
      timer.stop();
      double time = timer.time();

      // Output results of move statistics to files
      mcMoveManager_.output();
      if (Analyzer<D>::baseInterval > 0){
         analyzerManager_.output();
      }
      
      // Output how many times MDE has been solved for the simulation run
      Log::file() << std::endl;
      Log::file() << "MDE counter   " 
                  << system().compressor().counterMDE() << std::endl;
      Log::file() << std::endl;
      
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
   * Open, read and analyze a trajectory file
   */
   template <int D>
   void McSimulator<D>::analyzeTrajectory(int min, int max,
                                          std::string classname,
                                          std::string filename)
   {
      // Preconditions
      if (min < 0) UTIL_THROW("min < 0");
      if (max < 0) UTIL_THROW("max < 0");
      if (max < min) UTIL_THROW("max < min!");
      UTIL_CHECK(Analyzer<D>::baseInterval);
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
      Log::file() << "Reading " << filename << std::endl;
      trajectoryReaderPtr->open(filename);
      trajectoryReaderPtr->readHeader();
      // Read Header
      // Main loop over trajectory frames
      Timer timer;
      Log::file() << "Begin main loop" << std::endl;
      bool hasFrame = true;
      timer.start();
      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         hasFrame = trajectoryReaderPtr->readFrame();
         if (hasFrame) {
            clearData();
            // Initialize analyzers 
            if (iStep_ == min) analyzerManager_.setup();
            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager_.sample(iStep_);
               if ((iStep_ % 100) == 0){
                  Log::file() << "Analyzing steps: " << iStep_ << std::endl;
               }
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
      // Output time 
      Log::file() << std::endl;
      Log::file() << "# of frames   " << nFrames << std::endl;
      Log::file() << "run time      " << timer.time() 
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames) 
                  << "  sec" << std::endl;
      Log::file() << std::endl;
   }
   
   template<int D>
   void McSimulator<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "McSimulator times contributions:\n";
      mcMoveManager_.outputTimers(out);
   }
   
   template<int D>
   void McSimulator<D>::clearTimers()
   {
      mcMoveManager_.clearTimers();
   }

   /*
   * Save the current Monte-Carlo state.
   *
   * Used before attempting a Monte-Carlo move.
   */
   template <int D>
   void McSimulator<D>::saveMcState()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(hasWC_);
      UTIL_CHECK(hasHamiltonian_);
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(!mcState_.hasData);
      
      int nMonomer = system().mixture().nMonomer();
      int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      // Set field components
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>> (mcState_.w[i].cField(), 
             system().w().rgrid(i).cField(), meshSize);
         assignReal<<<nBlocks, nThreads>>>
            (mcState_.wc[i].cField(), wc_[i].cField(), meshSize);
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
   * Used when an attempted Monte-Carlo move is rejected.
   */
   template <int D>
   void McSimulator<D>::restoreMcState()
   {
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(mcState_.hasData);
      
      int nMonomer = system().mixture().nMonomer();
      int meshSize = system().domain().mesh().size();
      
      system().setWRGrid(mcState_.w); 

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (wc_[i].cField(), mcState_.wc[i].cField(), meshSize);
      }
      hamiltonian_ = mcState_.hamiltonian;
      idealHamiltonian_ = mcState_.idealHamiltonian;
      fieldHamiltonian_ = mcState_.fieldHamiltonian;
      hasHamiltonian_ = true;
      hasWC_ = true;
      mcState_.hasData = false;
   }
 
   /*
   * Clear the saved Monte-Carlo state.
   *
   * Used when an attempted Monte-Carlo move is accepted.
   */
   template <int D>
   void McSimulator<D>::clearMcState()
   {  mcState_.hasData = false; }

}
}
#endif
