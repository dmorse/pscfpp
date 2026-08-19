#ifndef RP_BD_SIMULATOR_H
#define RP_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>
#include <pscf/backends/TmplDeclare.h>

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
      template <int D, class T> class BdStep;
      template <int D, class T> class Analyzer;
      template <int D, class T> class AnalyzerManager;
      template <int D, class T> class TrajectoryReader;
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Brownian dynamics simulator for PS-FTS.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class CPT or CUT
   *
   * \see \ref rp_BdSimulator_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class BdSimulator : public Simulator<D,T>
   {

   public:

      /// \name Construction, destruction and initialization
      ///@{

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      BdSimulator(System<D,T>& system);

      /**
      * Destructor.
      */
      ~BdSimulator();

      // Prohibit copying or assignment.
      BdSimulator(BdSimulator<D,T> const &) = delete;
      BdSimulator<D,T>& operator = (BdSimulator<D,T> const &) = delete;

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
      BdStep<D,T>& bdStep();

      /**
      * Get the AnalyzerManager (const).
      */
      AnalyzerManager<D,T>& analyzerManager();

      /**
      * Get the AnalyzerManager (non-const).
      */
      AnalyzerManager<D,T> const& analyzerManager() const;

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory< TrajectoryReader<D,T> >& trajectoryReaderFactory();

      ///@}

   protected:

      /// Alias for Simulator class in program-level namespace.
      using SimulatorT = Simulator<D,T>;

      /// Alias for BdSimulator class in program-level namespace.
      using BdSimulatorT = BdSimulator<D,T>;

      // Inherited protected function (selected)
      using Simulator<D,T>::state;

      // Inherited protected data members (selected)
      using Simulator<D,T>::iStep_;
      using Simulator<D,T>::iTotalStep_;

   private:

      /// Private alias for analyzer class
      using AnalyzerT = Analyzer<D,T>;

      /**
      * Manager for Analyzer.
      */
      AnalyzerManager<D,T>* analyzerManagerPtr_;

      /**
      * Pointer to Brownian dynamics step algorithm.
      */
      BdStep<D,T>* bdStepPtr_;

      /**
      * Pointer to a BdStep factory.
      */
      Factory< BdStep<D,T> >* bdStepFactoryPtr_;

      /**
      * Pointer to a trajectory reader factory.
      */
      Factory< TrajectoryReader<D,T> >* trajectoryReaderFactoryPtr_;

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
   BdStep<D,T>& BdSimulator<D,T>::bdStep()
   {
      UTIL_CHECK(bdStepPtr_);
      return *bdStepPtr_;
   }

   // Get the AnalyzerManager (const).
   template <int D, class T> inline
   AnalyzerManager<D,T> const& 
   BdSimulator<D,T>::analyzerManager() const
   {  return *analyzerManagerPtr_; }

   // Get the AnalyzerManager (non-const).
   template <int D, class T> inline
   AnalyzerManager<D,T>& BdSimulator<D,T>::analyzerManager()
   {  return *analyzerManagerPtr_; }

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(BdSimulator)

}
}
#endif
