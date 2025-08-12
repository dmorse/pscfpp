#ifndef RPG_BD_SIMULATOR_TPP
#define RPG_BD_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"

#include <rpg/fts/brownian/BdStep.h>
#include <rpg/fts/brownian/BdStepFactory.h>
#include <rpg/fts/analyzer/AnalyzerFactory.h>
#include <rpg/fts/trajectory/TrajectoryReader.h>
#include <rpg/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpg/fts/perturbation/PerturbationFactory.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/fts/ramp/RampFactory.h>
#include <rpg/fts/ramp/Ramp.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>

#include <pscf/cuda/CudaRandom.h>
#include <util/misc/Timer.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BdSimulator<D>::BdSimulator(System<D>& system)
    : Simulator<D>(system),
      analyzerManager_(*this, system),
      bdStepPtr_(0),
      bdStepFactoryPtr_(0),
      trajectoryReaderFactoryPtr_(0),
      seed_(0)
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
   * Read instructions for creating objects from file.
   */
   template <int D>
   void BdSimulator<D>::readParameters(std::istream &in)
   {
      // Optionally read random seed, initialize random number generators
      readRandomSeed(in);

      // Optionally read and instantiate a BdStep object
      bool isEnd = false;
      std::string className;
      UTIL_CHECK(!hasBdStep());
      UTIL_CHECK(bdStepFactoryPtr_);
      bdStepPtr_ = bdStepFactoryPtr_->readObjectOptional(in, *this,
                                                         className, isEnd);

      // Attempt to read the optional Compressor, Perturbation and Ramp blocks
      readCompressor(in, isEnd);
      readPerturbation(in, isEnd);
      if (hasBdStep()) {
         UTIL_CHECK(hasCompressor());
         readRamp(in, isEnd);
      }
      // Compressor is required if a BdStep exists
      // A Ramp is allowed only if a BdStep exists

      // Read AnalyzerManager block
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);

      // Figure out what variables need to be saved in stored state
      state_.needsCc = false;
      state_.needsDc = false;
      state_.needsHamiltonian = false;
      if (hasBdStep()) {
         if (stepper().needsCc()){
            state_.needsCc = true;
         }
         if (stepper().needsDc()){
            state_.needsDc = true;
         }
      }

      // Initialize Simulator<D> base class
      allocate();
   }

   /*
   * Setup before the main loop of a simulate or analyze command.
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

      if (hasBdStep()) {
         stepper().setup();
      }

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
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(hasBdStep());
      UTIL_CHECK(hasCompressor());

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

      // Analysis initial step (if any)
      analyzerTimer.start();
      analyzerManager_.sample(iStep_);
      analyzerTimer.stop();

      for (iTotalStep_ = 0; iTotalStep_ < nStep; ++iTotalStep_) {

         // Take a step (modifies W fields)
         bool converged;
         converged = stepper().step();

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
            Log::file() << "Step: "<< iTotalStep_<< " fail to converge" << "\n";
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
         Log::file() << "nFail Step          " << (nStep - iStep_)
                     << std::endl;
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
   * Output timer results.
   */
   template<int D>
   void BdSimulator<D>::outputTimers(std::ostream& out) const
   {
      compressor().outputTimers(out);
   }

}
}
#endif
