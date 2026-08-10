#ifndef RP_BD_SIMULATOR_TPP
#define RP_BD_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdSimulator.h>

#include <rp/fts/brownian/BdStep.h>
#include <rp/fts/brownian/BdStepFactory.h>

#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/fts/analyzer/AnalyzerFactory.h>
#include <rp/fts/trajectory/TrajectoryReader.h>
#include <rp/fts/trajectory/TrajectoryReaderFactory.h>
#include <rp/fts/perturbation/Perturbation.h>
#include <rp/fts/ramp/Ramp.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <prdc/field/RField.h>

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
   BdSimulator<D,T>::BdSimulator(System<D,T>& system)
    : Simulator<D,T>(system),
      analyzerManagerPtr_(nullptr),
      bdStepPtr_(nullptr),
      bdStepFactoryPtr_(nullptr),
      trajectoryReaderFactoryPtr_(nullptr)
   {
      ParamComposite::setClassName("BdSimulator");
      analyzerManagerPtr_ 
             = new AnalyzerManager<D,T>(*this, system),
      bdStepFactoryPtr_ = new BdStepFactory<D,T>(*this);
      trajectoryReaderFactoryPtr_
             = new TrajectoryReaderFactory<D,T>(system);
      Analyzer<D,T>::initStatic();
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
      Simulator<D,T>::readRandomSeed(in);

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
      Simulator<D,T>::readCompressor(in, isEnd);
      if (hasBdStep()) {
         UTIL_CHECK((Simulator<D,T>::hasCompressor()));
      }

      // Optionally read a Perturbation block
      Simulator<D,T>::readPerturbation(in, isEnd);

      // Optionally read a Ramp block
      if (hasBdStep()) {
         Simulator<D,T>::readRamp(in, isEnd);
      }

      // Optionally read an AnalyzerManager block
      Analyzer<D,T>::baseInterval = 0; // default value
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

      // Allocate memory for Simulator<D,T> base class
      Simulator<D,T>::allocate();
   }

   /*
   * Setup before main loop of a simulate or analyze command.
   */
   template <int D, class T>
   void BdSimulator<D,T>::setup(int nStep)
   {
      UTIL_CHECK((Simulator<D,T>::system().w().hasData()));

      // Eigenanalysis of the projected chi matrix.
      Simulator<D,T>::analyzeChi();

      if (Simulator<D,T>::hasPerturbation()) {
         Simulator<D,T>::perturbation().setup();
      }

      if (Simulator<D,T>::hasRamp()) {
         Simulator<D,T>::ramp().setup(nStep);
      }

      // Solve MDE and compute c-fields for the intial state
      Simulator<D,T>::system().compute();

      // Compress the initial state (adjust pressure-like field)
      if (Simulator<D,T>::hasCompressor()) {
         Simulator<D,T>::compressor().compress();
         Simulator<D,T>::compressor().clearTimers();
      }

      // Compute field components and Hamiltonian for initial state.
      Simulator<D,T>::computeWc();
      Simulator<D,T>::computeCc();
      Simulator<D,T>::computeDc();
      Simulator<D,T>::computeHamiltonian();

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
      UTIL_CHECK((Simulator<D,T>::hasCompressor()));
      UTIL_CHECK((Simulator<D,T>::system().w().hasData()));

      // Initial setup
      setup(nStep);
      iStep_ = 0;
      if (Simulator<D,T>::hasRamp()) {
         Simulator<D,T>::ramp().setParameters(iStep_);
      }
      int analyzerBaseInterval = Analyzer<D,T>::baseInterval;

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

            if (Simulator<D,T>::hasRamp()) {
               Simulator<D,T>::ramp().setParameters(iStep_);
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
      if (Simulator<D,T>::hasRamp()){
         Log::file() << std::endl;
         Simulator<D,T>::ramp().output();
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
                  << Simulator<D,T>::compressor().mdeCounter() << std::endl;
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
      UTIL_CHECK((Analyzer<D,T>::baseInterval > 0));
      UTIL_CHECK(analyzerManager().size() > 0);

      // Construct TrajectoryReader
      TrajectoryReader<D,T>* trajectoryReaderPtr;
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
            Simulator<D,T>::clearData();

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
   Factory< TrajectoryReader<D,T> >& 
   BdSimulator<D,T>::trajectoryReaderFactory()
   {
      UTIL_CHECK(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

}
}
#endif
