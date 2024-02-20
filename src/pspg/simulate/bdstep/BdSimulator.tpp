#ifndef PSPG_BD_SIMULATOR_TPP
#define PSPG_BD_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"

#include <pspg/simulate/bdstep/BdStep.h>
#include <pspg/simulate/bdstep/BdStepFactory.h>
#include <pspg/simulate/analyzer/AnalyzerFactory.h>
#include <pspg/simulate/trajectory/TrajectoryReader.h>
#include <pspg/simulate/trajectory/TrajectoryReaderFactory.h>
#include <pspg/compressor/Compressor.h>
#include <pspg/System.h>

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
      // Read random seed. Default is seed_(0), taken from the clock time
      readOptional(in, "seed", seed_);
      // Set random seed
      random().setSeed(seed_);
      cudaRandom().setSeed(seed_);
      
      // Initialize Simulator<D> base class
      allocate();
      analyzeChi();

      // Instantiate an BdStep object
      std::string className;
      bool isEnd;
      bdStepPtr_ = 
          bdStepFactoryPtr_->readObject(in, *this, className, isEnd);

      // Read block of analyzer
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);

   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void BdSimulator<D>::setup()
   {
      UTIL_CHECK(system().w().hasData());
      analyzeChi();
      system().compute();
      computeWc();
      computeCc();
      computeDc();
      computeHamiltonian();
      stepper().setup();
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

      // Initial setup
      setup();

      // Main simulation loop
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

         // Take a step (modifies W fields)
         stepper().step();

      }

      // Analysis after final step (if any)
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

      // Output results analyzers to files
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
   * Output timer results.
   */
   template<int D>
   void BdSimulator<D>::outputTimers(std::ostream& out)
   {}

   /*
   * Clear all timers.
   */
   template<int D>
   void BdSimulator<D>::clearTimers()
   { }

}
}
#endif
