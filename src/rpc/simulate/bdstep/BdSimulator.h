#ifndef RPC_BD_SIMULATOR_H
#define RPC_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/simulate/Simulator.h>                  // base class
#include <rpc/simulate/analyzer/AnalyzerManager.h>   // member
#include <util/param/Factory.h>                      // member template

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class BdStep;
   template <int D> class TrajectoryReader;

   /**
   * Brownian dynamics simulator.
   *
   * \ingroup Rpc_Simulate_BdStep_Module
   */
   template <int D>
   class BdSimulator : public Simulator<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      BdSimulator(System<D>& system);

      /**
      * Destructor.
      */
      ~BdSimulator();

      /**
      * Read parameters for a MC simulation.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic Monte-Carlo simulation using the
      * partial saddle-point approximation.
      *
      * \param nStep  number of Monte-Carlo steps
      */
      void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader class
      * specified by the "classname" argument to read a trajectory
      * file.
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader class to use
      * \param filename  name of the trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Get BdStep.
      */
      BdStep<D>& stepper();

      /**
      * Get AnalyzerManger
      */
      AnalyzerManager<D>& analyzerManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory<TrajectoryReader<D>>& trajectoryReaderFactory();

      ///@}

      // Inherited public functions

      using Simulator<D>::system;
      using Simulator<D>::compressor;
      using Simulator<D>::random;
      using Simulator<D>::allocate;
      using Simulator<D>::analyzeChi;
      using Simulator<D>::computeWc;
      using Simulator<D>::computeCc;
      using Simulator<D>::computeDc;
      using Simulator<D>::wc;
      using Simulator<D>::cc;
      using Simulator<D>::dc;
      using Simulator<D>::hasWc;
      using Simulator<D>::hasCc;
      using Simulator<D>::hasDc;
      using Simulator<D>::clearData;
      using Simulator<D>::computeHamiltonian;
      using Simulator<D>::hamiltonian;
      using Simulator<D>::idealHamiltonian;
      using Simulator<D>::fieldHamiltonian;
      using Simulator<D>::hasHamiltonian;
      using Simulator<D>::hasPerturbation;
      using Simulator<D>::hasRamp;
      using Simulator<D>::saveState;
      using Simulator<D>::restoreState;
      using Simulator<D>::clearState;
      using Simulator<D>::outputMdeCounter;

   protected:

      // Inherited protected functions

      using ParamComponent::indent;
      using ParamComposite::setClassName;
      using ParamComposite::readParamComposite;
      using ParamComposite::readParamCompositeOptional;
      using ParamComposite::readOptional;
      using Simulator<D>::readCompressor;
      using Simulator<D>::perturbation;
      using Simulator<D>::perturbationFactory;
      using Simulator<D>::setPerturbation;
      using Simulator<D>::ramp;
      using Simulator<D>::rampFactory;
      using Simulator<D>::setRamp;

      // Inherited protected data members

      using Simulator<D>::wc_;
      using Simulator<D>::hasWc_;
      using Simulator<D>::hamiltonian_;
      using Simulator<D>::idealHamiltonian_;
      using Simulator<D>::fieldHamiltonian_;
      using Simulator<D>::hasHamiltonian_;
      using Simulator<D>::iStep_;
      using Simulator<D>::iTotalStep_;
      using Simulator<D>::state_;
      using Simulator<D>::seed_;

   private:

      /**
      * Manager for Analyzer.
      */
      AnalyzerManager<D> analyzerManager_;

      /**
      * Pointer to Brownian dynamics step algorithm.
      */
      BdStep<D>* bdStepPtr_;

      /**
      * Pointer to a trajectory reader/writer factory.
      */
      Factory< BdStep<D> >* bdStepFactoryPtr_;

      /**
      * Pointer to a trajectory reader/writer factory.
      */
      Factory< TrajectoryReader<D> >* trajectoryReaderFactoryPtr_;

      // Private member functions

      /**
      * Called at the beginning of the simulation.
      */
      void setup();

   };

   // Get the Brownian dynamics stepper.
   template <int D>
   inline BdStep<D>& BdSimulator<D>::stepper()
   {  return *bdStepPtr_; }

   // Get the analyzer manager.
   template <int D>
   inline AnalyzerManager<D>& BdSimulator<D>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReaderfactory
   template <int D>
   inline
   Factory<TrajectoryReader<D> >& BdSimulator<D>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

   #ifndef RPC_BD_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class BdSimulator<1>;
   extern template class BdSimulator<2>;
   extern template class BdSimulator<3>;
   #endif

}
}
#endif