#ifndef RP_BD_SIMULATOR_H
#define RP_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <iostream>
#include <string>

// Forward declaration
namespace Util {
   template <class T> class Factory;
}
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Brownian dynamics simulator for PS-FTS.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, both also named BdSimulator,
   * that are defined in Rpc and Rpg namespaces and used in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_BdSimulator_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class BdSimulator : public T::Simulator 
   {

   public:

      // Protected constructor and destructor (see below).

      // Not copyable or assignable.
      BdSimulator(BdSimulator<D,T> const &) = delete;
      BdSimulator<D,T>& operator = (BdSimulator<D,T> const &) = delete;

      /// \name Initialization
      ///@{

      /**
      * Read parameter file block.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      ///@}
      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Brownian dynamics (BD) simulation.
      *
      * Perform a field theoretic BD simulation using the partial
      * saddle-point approximation.
      *
      * \param nStep  number of BD steps
      */
      void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function creates an instance of the TrajectoryReader subclass
      * specified by the "classname" argument, and uses it to read and
      * analyze a section of a trajectory file, starting at frame number
      * "min" and ending at frame number "max".
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader subclass to use
      * \param filename  name of the trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Does this BdSimulator have a BdStep object?
      */
      bool hasBdStep() const;

      /**
      * Get the BdStep by non-const reference.
      */
      typename T::BdStep& bdStep();

      /**
      * Get the AnalyzerManager (const).
      */
      typename T::AnalyzerManager& analyzerManager();

      /**
      * Get the AnalyzerManager (non-const).
      */
      typename T::AnalyzerManager const& analyzerManager() const;

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory< typename T::TrajectoryReader >& trajectoryReaderFactory();

      ///@}

   protected:

      /// Alias for System class in program-level namespace.
      using SystemT = System<D,T>;

      /// Alias for Simulator class in program-level namespace.
      using SimulatorT = typename T::Simulator;

      /// Alias for BdSimulator class in program-level namespace.
      using BdSimulatorT = typename T::BdSimulator;

      /**
      * Constructor.
      *
      * \param system  parent System
      * \param bdSimulator  instance of enclosing subclass
      */
      BdSimulator(SystemT& system, BdSimulatorT& bdSimulator);

      /**
      * Destructor.
      */
      ~BdSimulator();

      // Inherited protected function (selected)
      using SimulatorT::state;

      // Inherited protected data members (selected)
      using SimulatorT::iStep_;
      using SimulatorT::iTotalStep_;

   private:

      /// Private alias for analyzer class
      using AnalyzerT = typename T::Analyzer;

      /**
      * Manager for Analyzer.
      */
      typename T::AnalyzerManager* analyzerManagerPtr_;

      /**
      * Pointer to Brownian dynamics step algorithm.
      */
      typename T::BdStep* bdStepPtr_;

      /**
      * Pointer to a BdStep factory.
      */
      Factory< typename T::BdStep >* bdStepFactoryPtr_;

      /**
      * Pointer to a trajectory reader factory.
      */
      Factory< typename T::TrajectoryReader >* trajectoryReaderFactoryPtr_;

      // Private member function

      /**
      * Setup before the main loop.
      *
      * \param nStep  number of steps planned for the simulation
      */
      void setup(int nStep);

   };

   // Inline member functions

   // Does this BdSimulator have a BdStep?
   template <int D, class T> inline
   bool BdSimulator<D,T>::hasBdStep() const
   {  return (bool)bdStepPtr_; }

   // Get the BdStep.
   template <int D, class T> inline
   typename T::BdStep& BdSimulator<D,T>::bdStep()
   {
      UTIL_CHECK(bdStepPtr_);
      return *bdStepPtr_;
   }

   // Get the AnalyzerManager (const).
   template <int D, class T> inline
   typename T::AnalyzerManager const& 
   BdSimulator<D,T>::analyzerManager() const
   {  return *analyzerManagerPtr_; }

   // Get the AnalyzerManager (non-const).
   template <int D, class T> inline
   typename T::AnalyzerManager& BdSimulator<D,T>::analyzerManager()
   {  return *analyzerManagerPtr_; }

}
}
#endif
