#ifndef PSPC_MC_SIMULATOR_H
#define PSPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/simulate/Simulator.h>                   // member
#include <pspc/simulate/mcmove/McState.h>              // member
#include <pspc/simulate/mcmove/McMoveManager.h>        // member
#include <pspc/simulate/analyzer/AnalyzerManager.h>    // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   template <int D> class McMove;
   template <int D> class TrajectoryReader;

   /**
   * Monte-Carlo simulation coordinator.
   *
   * \ingroup Rpc_Simulate_Module
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

      ///@}
      /// \name Utilities for MC Moves
      ///@{

      /**
      * Save a copy of the Monte-Carlo state.
      *
      * This function and restoreMcState() are intended for use in
      * the implementation of field theoretic Monte Carlo moves. This
      * function stores the current w fields and the corresponding
      * Hamiltonian value.  This is normally the first step of a MC
      * move, prior to an attempted modification of the fields stored
      * in the system w field container.
      */
      void saveMcState();

      /**
      * Restore the saved copy of the Monte-Carlo state.
      *
      * This function  and saveMcState() are intended to be used
      * together in the implementation of Monte-Carlo moves. If an
      * attempted move is rejected, restoreMcState() is called to
      * restore the fields ahd Hamiltonian value that were saved
      * by a previous call to the function saveMcState().
      */
      void restoreMcState();

      /**
      * Clear the saved copy of the Monte-Carlo state.
      *
      * This function, restoreMcState(), and saveMcState() are intended
      * to be used together in the implementation of Monte-Carlo moves. If
      * an attempted move is accepted, clearMcState() is called to clear
      * clear mcState_.hasData
      */
      void clearMcState();

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
      using Simulator<D>::wc;
      using Simulator<D>::hasWc;
      using Simulator<D>::clearData;
      using Simulator<D>::computeHamiltonian;
      using Simulator<D>::hamiltonian;
      using Simulator<D>::idealHamiltonian;
      using Simulator<D>::fieldHamiltonian;
      using Simulator<D>::hasHamiltonian;

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
      * State saved during MC simulation.
      */
      mutable McState<D> mcState_;

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

   #ifndef PSPC_MC_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class McSimulator<1>;
   extern template class McSimulator<2>;
   extern template class McSimulator<3>;
   #endif

}
}
#endif
