#ifndef RPC_MC_SIMULATOR_TPP
#define RPC_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <rpc/system/System.h>
#include <rpc/fts/montecarlo/McMoveFactory.h>
#include <rpc/fts/analyzer/AnalyzerFactory.h>
#include <rpc/fts/trajectory/TrajectoryReader.h>
#include <rpc/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/fts/perturbation/PerturbationFactory.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/ramp/RampFactory.h>
#include <rpc/fts/ramp/Ramp.h>

#include <util/random/Random.h>
#include <util/misc/Timer.h>
#include <util/global.h>

#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McSimulator<D>::McSimulator(System<D>& system)
    : Simulator<D>(system),
      mcMoveManager_(*this, system),
      analyzerManager_(*this, system),
      trajectoryReaderFactoryPtr_(nullptr)
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
   * Read parameter file block.
   */
   template <int D>
   void McSimulator<D>::readParameters(std::istream &in)
   {
      // Read optional random seed value
      readRandomSeed(in);

      // Read optional McMoveManager block
      readParamCompositeOptional(in, mcMoveManager_);

      bool isEnd = false;

      // Read optionally Compressor block
      readCompressor(in, isEnd);
      if (hasMcMoves()) {
         UTIL_CHECK(hasCompressor());
      }

      // Read optional Perturbation and/or Ramp blocks
      readPerturbation(in, isEnd);
      if (hasMcMoves()) {
         readRamp(in, isEnd);
      }

      // Read optional AnalyzerManager block
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);

      // Figure out what needs to be saved
      state_.needsCc = false;
      state_.needsDc = false;
      state_.needsHamiltonian = true;
      if (mcMoveManager_.needsCc()){
         state_.needsCc = true;
      }
      if (mcMoveManager_.needsDc()){
         state_.needsDc = true;
      }

      // Initialize Simulator<D> base class
      Simulator<D>::allocate();

   }

   /*
   * Initialize just prior to a run.
   */
   template <int D>
   void McSimulator<D>::setup(int nStep)
   {
      UTIL_CHECK(system().w().hasData());

      // Eigenanalysis of the projected chi matrix.
      analyzeChi();
      
      if (hasPerturbation()) {
         perturbation().setup();
      }
      
      if (hasRamp()) {
         ramp().setup(nStep);
      }
   
      // Solve MDE and compute c-fields for the intial state
      system().compute();
      
      // Compress the initial state (adjust pressure-like field)
      if (hasCompressor()) {
         compressor().compress();
         compressor().clearTimers();
      }
      
      // Compute field components and Hamiltonian
      computeWc();
      if (state_.needsCc || state_.needsDc) {
         computeCc();
      }
      if (state_.needsDc) {
         computeDc();
      }
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
      UTIL_CHECK(hasMcMoves());
      UTIL_CHECK(hasCompressor());
      
      // Initial setup
      setup(nStep);
   
      Log::file() << std::endl;

      // Initialize timer, iStep_, and ramp (if any)
      Timer timer;
      Timer analyzerTimer;
      timer.start();
      iStep_ = 0;
      if (hasRamp()) {
         ramp().setParameters(iStep_);
      }

      // Analysis initial step (if any)
      analyzerTimer.start();
      analyzerManager_.sample(iStep_);
      analyzerTimer.stop();

      // Main Monte Carlo loop
      for (iTotalStep_ = 0; iTotalStep_ < nStep; ++iTotalStep_) {

         // Choose and attempt an McMove
         bool converged;
         converged = mcMoveManager_.chooseMove().move();

         if (converged){
            iStep_++;
            
            if (hasRamp()) {
               ramp().setParameters(iStep_);               
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

         } else{
            Log::file() << "Step: "<< iTotalStep_ 
                        << " failed to converge" << "\n";
         }
      }

      timer.stop();
      double time = timer.time();
      double analyzerTime = analyzerTimer.time();

      // Output results of move statistics to files
      mcMoveManager_.output();
      if (Analyzer<D>::baseInterval > 0){
         analyzerManager_.output();
      }
      
      // Output results of ramp
      if (hasRamp()){
         ramp().output();
      }

      // Output times for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep               " << nStep << std::endl;
      if (iStep_ != nStep){
         Log::file() << "nFail Step          " << (nStep - iStep_) << std::endl;
      }
      Log::file() << "Total run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep        " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << "Analyzer run time   " << analyzerTime
                  << " sec" << std::endl;
      Log::file() << std::endl;

      // Output number of times MDE has been solved for the simulation run
      outputMdeCounter(Log::file());
      
      // Print McMove acceptance statistics
      long attempt;
      long accept;
      long fail;
      using namespace std;
      Log::file() << "Move Statistics:" << endl << endl;
      Log::file() << setw(20) << left <<  "Move Name"
           << setw(10) << right << "Attempted"
           << setw(10) << right << "Accepted"
           << setw(13) << right << "AcceptRate"
           << setw(10) << right << "Failed"
           << setw(13) << right << "FailRate"
           << endl;
      int nMove = mcMoveManager_.size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = mcMoveManager_[iMove].nAttempt();
         accept  = mcMoveManager_[iMove].nAccept();
         fail  = mcMoveManager_[iMove].nFail();
         Log::file() << setw(20) << left
              << mcMoveManager_[iMove].className()
              << setw(10) << right << attempt
              << setw(10) << accept
              << setw(13) << fixed << setprecision(5)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << setw(10) << fail
              << setw(13) << fixed << setprecision(5)
              << ( attempt == 0 ? 0.0 : double(fail)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

   }

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
      bool hasFrame;
      timer.start();
      hasFrame = trajectoryReaderPtr->readFrame();
      
      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         if (hasFrame) {
            clearData();

            // Initialize analyzers
            if (iStep_ == min) {
               setup(iStep_);
            }

            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager_.sample(iStep_);
            }
         }
         
         hasFrame = trajectoryReaderPtr->readFrame();
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
   void McSimulator<D>::outputTimers(std::ostream& out) const
   {
      compressor().outputTimers(out);
      out << "\n";
      out << "MC move time contributions:\n";
      mcMoveManager_.outputTimers(out);
      out << "\n";
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
