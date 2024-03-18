#ifndef RPC_MC_SIMULATOR_H
#define RPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/simulate/Simulator.h>                   // member
#include <rpc/simulate/mcmove/McMoveManager.h>        // member
#include <rpc/simulate/analyzer/AnalyzerManager.h>    // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   template <int D> class McMove;
   template <int D> class TrajectoryReader;

   /**
   * Monte-Carlo simulation coordinator.
   *
   * \ingroup Rpc_Simulate_McMove_Module
   */
   template <int D>
   class McSimulator : public Simulator<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McSimulator(System<D>& system);

      /**
      * Destructor.
      */
      ~McSimulator();

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

      /**
      * Output timing results
      */
      virtual void outputTimers(std::ostream& out);

      /**
      * Clear timers
      */
      virtual void clearTimers();
      
      /**
      * Return the simulations whether needs to store cc fields
      */
      bool needsCc();
      
      /**
      * Return the simulations whether needs to store Dc fields
      */
      bool needsDc();

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Get McMoveManger
      */
      McMoveManager<D>& mcMoveManager();

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
      using Simulator<D>::random;
      using Simulator<D>::allocate;
      using Simulator<D>::analyzeChi;
      using Simulator<D>::computeWc;
      using Simulator<D>::computeCc;
      using Simulator<D>::computeDc;
      using Simulator<D>::wc;
      using Simulator<D>::hasWc;
      using Simulator<D>::clearData;
      using Simulator<D>::computeHamiltonian;
      using Simulator<D>::hamiltonian;
      using Simulator<D>::idealHamiltonian;
      using Simulator<D>::fieldHamiltonian;
      using Simulator<D>::hasHamiltonian;
      using Simulator<D>::saveState;
      using Simulator<D>::restoreState;
      using Simulator<D>::clearState;

   protected:

      // Inherited protected functions

      using ParamComposite::setClassName;
      using ParamComposite::readParamComposite;
      using ParamComposite::readParamCompositeOptional;
      using ParamComposite::readOptional;

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

   private:

      /**
      * Manger for Monte Carlo Move.
      */
      McMoveManager<D> mcMoveManager_;

      /**
      * Manger for Monte Carlo Analyzer.
      */
      AnalyzerManager<D> analyzerManager_;
      
      /**
      * Pointer to a trajectory reader/writer factory.
      */
      Factory< TrajectoryReader<D> >* trajectoryReaderFactoryPtr_;

      /**
      * Seed value of random generator. 
      */
      long seed_;
      
      // Private member functions

      /**
      * Called at the beginning of the simulation member function.
      */
      void setup();

   };

   // Get the Monte-Carlo move manager.
   template <int D>
   inline McMoveManager<D>& McSimulator<D>::mcMoveManager()
   {  return mcMoveManager_; }

   // Get the Monte-Carlo analyzer manager.
   template <int D>
   inline AnalyzerManager<D>& McSimulator<D>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReaderfactory
   template <int D>
   inline 
   Factory<TrajectoryReader<D> >& McSimulator<D>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }
   
   // Return the simulations whether needs to store Cc fields
   template <int D>
   inline bool McSimulator<D>::needsCc()
   {
      return mcMoveManager_.needsCc();
   }
   
   // Return the simulations whether needs to store Dc fields
   template <int D>
   inline bool McSimulator<D>::needsDc()
   {
      return mcMoveManager_.needsDc();
   }
   
   #ifndef RPC_MC_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class McSimulator<1>;
   extern template class McSimulator<2>;
   extern template class McSimulator<3>;
   #endif

}
}
#endif
