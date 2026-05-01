#ifndef RP_BD_SIMULATOR_TPP
#define RP_BD_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"

#include <util/param/Factory.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/misc/Timer.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   BdSimulator<D,T>::BdSimulator(SystemT& system, BdSimulatorT& bdSimulator)
    : SimulatorT(system),
      analyzerManagerPtr_(nullptr),
      bdStepPtr_(nullptr),
      bdStepFactoryPtr_(nullptr),
      trajectoryReaderFactoryPtr_(nullptr)
   {
      ParamComposite::setClassName("BdSimulator");
      analyzerManagerPtr_ 
             = new typename T::AnalyzerManager(bdSimulator, system),
      bdStepFactoryPtr_ = new typename T::BdStepFactory(bdSimulator);
      trajectoryReaderFactoryPtr_
             = new typename T::TrajectoryReaderFactory(system);
      AnalyzerT::initStatic();
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   BdSimulator<D,T>::~BdSimulator()
   {
      delete analyzerManagerPtr_;
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
   template <int D, class T>
   void BdSimulator<D,T>::readParameters(std::istream &in)
   {
      // Optionally read random seed, initialize random number generators
      SimulatorT::readRandomSeed(in);

      // Optionally read a BdStep block
      bool isEnd = false;
      std::string className;
      UTIL_CHECK(!hasBdStep());
      UTIL_CHECK(bdStepFactoryPtr_);
      bdStepPtr_ =
         bdStepFactoryPtr_->readObjectOptional(in, *this,
                                               className, isEnd);
      if (!hasBdStep() && ParamComponent::echo()) {
         Log::file() << ParamComponent::indent() 
                     << "  BdStep{ [absent] }\n";
      }

      // Compressor is required if a BdStep exists
      // A Ramp is allowed only if a BdStep exists

      // Optionally read a Compressor block
      SimulatorT::readCompressor(in, isEnd);
      if (hasBdStep()) {
         UTIL_CHECK(SimulatorT::hasCompressor());
      }

      // Optionally read a Perturbation block
      SimulatorT::readPerturbation(in, isEnd);

      // Optionally read a Ramp block
      if (hasBdStep()) {
         SimulatorT::readRamp(in, isEnd);
      }

      // Optionally read an AnalyzerManager block
      AnalyzerT::baseInterval = 0; // default value
      ParamComposite::readParamCompositeOptional(in, analyzerManager());

      // Figure out what variables need to be saved in stored state()
      state().needsCc = false;
      state().needsDc = false;
      state().needsHamiltonian = false;
      if (hasBdStep()) {
         if (bdStep().needsCc()){
            state().needsCc = true;
         }
         if (bdStep().needsDc()){
            state().needsDc = true;
         }
      }

      // Allocate memory for SimulatorT base class
      SimulatorT::allocate();
   }

   /*
   * Setup before main loop of a simulate or analyze command.
   */
   template <int D, class T>
   void BdSimulator<D,T>::setup(int nStep)
   {
      UTIL_CHECK(SimulatorT::system().w().hasData());

      // Eigenanalysis of the projected chi matrix.
      SimulatorT::analyzeChi();

      if (SimulatorT::hasPerturbation()) {
         SimulatorT::perturbation().setup();
      }

      if (SimulatorT::hasRamp()) {
         SimulatorT::ramp().setup(nStep);
      }

      // Solve MDE and compute c-fields for the intial state
      SimulatorT::system().compute();

      // Compress the initial state (adjust pressure-like field)
      if (SimulatorT::hasCompressor()) {
         SimulatorT::compressor().compress();
         SimulatorT::compressor().clearTimers();
      }

      // Compute field components and Hamiltonian for initial state.
      SimulatorT::computeWc();
      SimulatorT::computeCc();
      SimulatorT::computeDc();
      SimulatorT::computeHamiltonian();

      if (hasBdStep()) {
         bdStep().setup();
      }

      if (analyzerManager().size() > 0){
         analyzerManager().setup();
      }

   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D, class T>
   void BdSimulator<D,T>::simulate(int nStep)
   {
      UTIL_CHECK(hasBdStep());
      UTIL_CHECK(SimulatorT::hasCompressor());
      UTIL_CHECK(SimulatorT::system().w().hasData());

      // Initial setup
      setup(nStep);
      iStep_ = 0;
      if (SimulatorT::hasRamp()) {
         SimulatorT::ramp().setParameters(iStep_);
      }
      int analyzerBaseInterval = AnalyzerT::baseInterval;

      // Start timer
      Timer timer;
      Timer analyzerTimer;
      timer.start();

      // Analysis for initial state (if any)
      analyzerTimer.start();
      if (analyzerManager().size() > 0) {
         analyzerManager().sample(iStep_);
      }
      analyzerTimer.stop();

      // Main simulation loop
      for (iTotalStep_ = 0; iTotalStep_ < nStep; ++iTotalStep_) {

         // Take a step (modifies W fields, then applies compressor)
         bool converged;
         converged = bdStep().step();

         // Accept step iff compressor converged
         if (converged){
            iStep_++;

            if (SimulatorT::hasRamp()) {
               SimulatorT::ramp().setParameters(iStep_);
            }

            // Analysis (if any)
            analyzerTimer.start();
            if (analyzerBaseInterval != 0) {
               if (analyzerManager().size() > 0) {
                  if (iStep_ % analyzerBaseInterval == 0) {
                     analyzerManager().sample(iStep_);
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

      // Output final analyzer results
      if (analyzerBaseInterval != 0){
         analyzerManager().output();
      }

      // Output results of ramp
      if (SimulatorT::hasRamp()){
         Log::file() << std::endl;
         SimulatorT::ramp().output();
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
                  << SimulatorT::compressor().mdeCounter() << std::endl;
      Log::file() << std::endl;

   }

   /*
   * Open, read and analyze a trajectory file
   */
   template <int D, class T>
   void BdSimulator<D,T>::analyze(int min,
                                  int max,
                                  std::string classname,
                                  std::string filename)
   {
      // Preconditions
      UTIL_CHECK(min >= 0);
      UTIL_CHECK(max >= min);
      UTIL_CHECK(AnalyzerT::baseInterval > 0);
      UTIL_CHECK(analyzerManager().size() > 0);

      // Construct TrajectoryReader
      typename T::TrajectoryReader* trajectoryReaderPtr;
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
      timer.start();
      bool hasFrame = trajectoryReaderPtr->readFrame();

      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         if (hasFrame) {
            SimulatorT::clearData();

            // Initialize analyzers
            if (iStep_ == min) {
               setup(iStep_);
            }

            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager().sample(iStep_);
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
      analyzerManager().output();

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
   * Get the TrajectoryReader factory.
   */
   template <int D, class T>
   Factory< typename T::TrajectoryReader >& 
   BdSimulator<D,T>::trajectoryReaderFactory()
   {
      UTIL_CHECK(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

}
}
#endif
