#ifndef RPC_BD_SIMULATOR_TPP
#define RPC_BD_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"

#include <rpc/fts/brownian/BdStep.h>
#include <rpc/fts/brownian/BdStepFactory.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/fts/analyzer/AnalyzerFactory.h>
#include <rpc/fts/trajectory/TrajectoryReader.h>
#include <rpc/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpc/fts/perturbation/PerturbationFactory.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/ramp/RampFactory.h>
#include <rpc/fts/ramp/Ramp.h>
#include <rpc/System.h>

#include <util/random/Random.h>
#include <util/misc/Timer.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BdSimulator<D>::BdSimulator(System<D>& system)
    : Simulator<D>(system),
      analyzerManager_(*this, system),
      bdStepPtr_(nullptr),
      bdStepFactoryPtr_(nullptr),
      trajectoryReaderFactoryPtr_(nullptr)
   {
      setClassName("BdSimulator");
      bdStepFactoryPtr_ = new BdStepFactory<D>(*this);
      trajectoryReaderFactoryPtr_
             = new TrajectoryReaderFactory<D>(system);
   }

   /*
   * Destructor.
   */
   template <int D>
   BdSimulator<D>::~BdSimulator()
   {
      if (bdStepFactoryPtr_) {
         delete bdStepFactoryPtr_;
      }
      if (bdStepPtr_) {
         delete bdStepPtr_;
      }
      if (trajectoryReaderFactoryPtr_) {
         delete trajectoryReaderFactoryPtr_;
      }
   }

   /*
   * Read parameter file block for a BD simulator.
   */
   template <int D>
   void BdSimulator<D>::readParameters(std::istream &in)
   {
      // Optionally read a random seed value
      readRandomSeed(in);
   
      // Optionally read a BdStep 
      bool isEnd = false;
      std::string className;
      bdStepPtr_ =
         bdStepFactoryPtr_->readObjectOptional(in, *this, 
                                               className, 
                                               isEnd);
      if (!hasBdStep() && ParamComponent::echo()) {
         Log::file() << indent() << "  BdStep{ [absent] }\n";
      }

      // Optionally read Compressor, Perturbation and/or Ramp blocks
      readCompressor(in, isEnd);
      if (hasBdStep()) {
         UTIL_CHECK(hasCompressor());
      }

      // Optionally Perturbation and/or Ramp blocks
      readPerturbation(in, isEnd);
      if (hasBdStep()) {
         readRamp(in, isEnd);
      }

      // Optionally read an AnalyzerManager
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);

      // Figure out what variables need to be saved
      state_.needsCc = false;
      state_.needsDc = false;
      state_.needsHamiltonian = false;
      if (hasBdStep()) {
         if (bdStep().needsCc()){
            state_.needsCc = true;
         }
         if (bdStep().needsDc()){
            state_.needsDc = true;
         }
      }

      // Allocate memory for Simulator<D> base class
      Simulator<D>::allocate();

   }

   /*
   * Setup before main loop of a BD simulation.
   */
   template <int D>
   void BdSimulator<D>::setup(int nStep)
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

      // Compute field components and Hamiltonian for initial state.
      computeWc();
      computeCc();
      computeDc();
      computeHamiltonian();

      bdStep().setup();
      if (analyzerManager_.size() > 0){
         analyzerManager_.setup();
      }

   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void BdSimulator<D>::simulate(int nStep)
   {
      UTIL_CHECK(hasBdStep());
      UTIL_CHECK(hasCompressor());
      UTIL_CHECK(system().w().hasData());

      // Initial setup
      setup(nStep);

      // Main simulation loop
      Timer timer;
      Timer analyzerTimer;
      timer.start();
      iStep_ = 0;
      if (hasRamp()) {
         ramp().setParameters(iStep_);
      }

      // Analysis for initial state (if any)
      analyzerTimer.start();
      if (analyzerManager_.size() > 0){
         analyzerManager_.sample(iStep_);
      }
      analyzerTimer.stop();

      for (iTotalStep_ = 0; iTotalStep_ < nStep; ++iTotalStep_) {

         // Take a step (modifies W fields)
         bool converged;
         converged = bdStep().step();

         if (converged){
            iStep_++;

            if (hasRamp()) {
               ramp().setParameters(iStep_);
            }

            // Analysis (if any)
            analyzerTimer.start();
            if (Analyzer<D>::baseInterval != 0) {
               if (analyzerManager_.size() > 0) {
                  if (iStep_ % Analyzer<D>::baseInterval == 0) {
                     analyzerManager_.sample(iStep_);
                  }
               }
            }
            analyzerTimer.stop();

         } else {
            Log::file() << "Step: "<< iTotalStep_
                        << " failed to converge" << "\n";
         }

      }

      timer.stop();
      double time = timer.time();
      double analyzerTime = analyzerTimer.time();

      // Output results analyzers to files
      if (Analyzer<D>::baseInterval > 0){
         analyzerManager_.output();
      }

      // Output results of ramp
      if (hasRamp()){
         Log::file() << std::endl;
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
      Log::file() << "MDE counter   "
                  << compressor().mdeCounter() << std::endl;
      Log::file() << std::endl;

      // Output compressor timer results
      // compressor().outputTimers(Log::file());
   }

   /*
   * Open, read and analyze a trajectory file
   */
   template <int D>
   void BdSimulator<D>::analyze(int min, int max,
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
               //analyzerManager_.setup();
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

}
}
#endif
