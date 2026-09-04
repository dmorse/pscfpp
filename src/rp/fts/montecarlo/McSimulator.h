#ifndef RP_MC_SIMULATOR_H
#define RP_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h> // base class

#include <pscf/backend/TmplDeclare.h>
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
      template <int D, class T> class McMoveManager; 
      template <int D, class T> class AnalyzerManager; 
      template <int D, class T> class TrajectoryReader; 
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Monte-Carlo simulation coordinator.
   * 
   * An McSimulator provides tools that are specific to MC simulation, 
   * in addition to those inherited from the Simulator base class. An
   * McSimulator has an McMoveManager and an AnalyzerManager, in addition
   * data members inherited from the Simulator class.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named McSimulator, that
   * are defined in Rpc and Rpg namespaces for use in the pscf_rpc and
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_McSimulator_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class McSimulator : public Simulator<D,T>
   {

   public:

      /// \name Construction, destruction and initialization
      ///@{

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      McSimulator(System<D,T>& system);

      /**
      * Destructor.
      */
      ~McSimulator();

      // Prohibit copying and assignment. 
      McSimulator(McSimulator<D,T> const &) = delete;
      McSimulator<D,T>& operator = (McSimulator<D,T> const &) = delete;

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
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic Monte-Carlo simulation using the
      * partial saddle-point approximation.
      *
      * \param nStep  number of attempted Monte-Carlo steps
      */
      void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader subclass
      * specified by the "classname" argument to read a trajectory
      * file.
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader subclass to use
      * \param filename  name of the trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      /**
      * Output timing results
      */
      virtual void outputTimers(std::ostream& out) const;

      /**
      * Clear timers
      */
      virtual void clearTimers();

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Get the McMoveManager (const).
      */
      McMoveManager<D,T> const& mcMoveManager() const;

      /**
      * Get the McMoveManager (non-const).
      */
      McMoveManager<D,T>& mcMoveManager();

      /**
      * Get the AnalyzerManager (const).
      */
      AnalyzerManager<D,T> const & analyzerManager() const;

      /**
      * Get the AnalyzerManager (non-const).
      */
      AnalyzerManager<D,T>& analyzerManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory<TrajectoryReader<D,T> >& trajectoryReaderFactory();

      /**
      * Have any MC moves been defined?
      *
      * Equivalent to a test for mcMoveManager().size() > 0.
      */
      bool hasMcMoves() const;

      /**
      * Does the stored state need to include Cc fields?
      *
      * This returns if McMove::needsCc returns true for one or more
      * MC moves ih the McMoveManager. Most MC moves require storage
      * of components of the c fields.
      */
      bool needsCc();

      /**
      * Does the stored state need to include Dc fields?
      *
      * This returns if McMove::needsDc returns true for one or more MC
      * moves in the McMoveManager. Only moves that uses forces to 
      * generate a proposed change in the fields, such as a "smart" MC 
      * force bias move, will require storage of these fields.
      */
      bool needsDc();

      ///@}

   protected:

      // Inherited protected member function (selected).
      using Simulator<D,T>::state;

      // Inherited protected data members (selected).
      using Simulator<D,T>::iStep_;
      using Simulator<D,T>::iTotalStep_;

   private:

      /**
      * Pointer to Manager for Monte Carlo moves.
      */
      McMoveManager<D,T>* mcMoveManagerPtr_;

      /**
      * Pointer to Manager for analyzers.
      */
      AnalyzerManager<D,T>* analyzerManagerPtr_;

      /**
      * Pointer to a trajectory reader factory.
      */
      Factory<TrajectoryReader<D,T> >* trajectoryReaderFactoryPtr_;

      // Private member function

      /**
      * Setup before the main loop of a simulate or analyze command.
      *
      * \param nStep  number of MC steps to attempt
      */
      void setup(int nStep);

   };

   // Inline member function definitions

   // Get the Monte-Carlo move manager (const)
   template <int D, class T> inline
   McMoveManager<D,T> const& McSimulator<D,T>::mcMoveManager() const
   {  
      UTIL_ASSERT(mcMoveManagerPtr_);
      return *mcMoveManagerPtr_; 
   }

   // Get the Monte-Carlo move manager (non-const).
   template <int D, class T> inline
   McMoveManager<D,T>& McSimulator<D,T>::mcMoveManager()
   {  
      UTIL_ASSERT(mcMoveManagerPtr_);
      return *mcMoveManagerPtr_; 
   }

   // Get the Analyzer manager (const).
   template <int D, class T> inline
   AnalyzerManager<D,T> const& McSimulator<D,T>::analyzerManager() 
   const
   {  
      UTIL_ASSERT(analyzerManagerPtr_);
      return *analyzerManagerPtr_; 
   }

   // Get the Analyzer manager.
   template <int D, class T> inline
   AnalyzerManager<D,T>& McSimulator<D,T>::analyzerManager()
   {  
      UTIL_ASSERT(analyzerManagerPtr_);
      return *analyzerManagerPtr_; 
   }

   // Get the TrajectoryReader factory.
   template <int D, class T> inline
   Factory<TrajectoryReader<D,T> >& 
   McSimulator<D,T>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

   // Have any MC moves been defined?
   template <int D, class T> inline
   bool McSimulator<D,T>::hasMcMoves() const
   {  return (bool)(mcMoveManager().size() > 0); }

   // Does the stored state need to include Cc fields?
   template <int D, class T> inline
   bool McSimulator<D,T>::needsCc()
   {  return state().needsCc; }

   // Does the stored state need to include Dc fields?
   template <int D, class T> inline
   bool McSimulator<D,T>::needsDc()
   {  return state().needsDc; }

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(McSimulator)

}
}
#endif
