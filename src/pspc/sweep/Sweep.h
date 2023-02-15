#ifndef PSPC_SWEEP_H
#define PSPC_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/sweep/SweepTmpl.h>          // base class template
#include <pspc/sweep/BasisFieldState.h>    // base class template parameter
#include "SweepParameter.h" // parameter class
#include <util/global.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Solve a sequence of problems along a line in parameter space.
   */
   template <int D>
   class Sweep : public SweepTmpl< BasisFieldState<D> >
   {

   public:

      /**
      * Default Constructor.
      */
      Sweep();

      /**
      * Constructor, creates assocation with parent system.
      */
      Sweep(System<D>& system);

      /**
      * Destructor.
      */
      ~Sweep();

      /**
      * Set association with parent System.
      */
      void setSystem(System<D>& system);

      /**
      * Read parameters from param file.
      * 
      * \param in Input stream from param file.
      */
      virtual void readParameters(std::istream& in);

      // Public members inherited from base class template SweepTmpl
      using SweepTmpl< BasisFieldState<D> >::historyCapacity;
      using SweepTmpl< BasisFieldState<D> >::historySize;
      using SweepTmpl< BasisFieldState<D> >::nAccept;
      using SweepTmpl< BasisFieldState<D> >::state;
      using SweepTmpl< BasisFieldState<D> >::s;
      using SweepTmpl< BasisFieldState<D> >::c;

   protected:

      /**
      * Check allocation state of fields in one state, allocate if necessary.
      *
      * \param state object that represents a stored state of the system.
      */
      virtual void checkAllocation(BasisFieldState<D>& state);

      /**
      * Setup operation at the beginning of a sweep.
      */
      virtual void setup();

      /**
      * Set non-adjustable system parameters to new values.
      *
      * \param sNew contour variable value for new trial solution.
      */
      virtual void setParameters(double sNew) = 0;

      /**
      * Create a guess for adjustable variables by continuation.
      *
      * \param sNew contour variable value for new trial solution.
      */
      virtual void extrapolate(double sNew);

      /**
      * Call current iterator to solve SCFT problem.
      *
      * Return 0 for sucessful solution, 1 on failure to converge.
      */
      virtual int solve(bool isContinuation);

      /**
      * Reset system to previous solution after iterature failure.
      *
      * The implementation of this function should reset the system state
      * to correspond to that stored in state(0).
      */
      virtual void reset();

      /**
      * Update state(0) and output data after successful convergence
      *
      * The implementation of this function should copy the current 
      * system state into state(0) and output any desired information
      * about the current converged solution.
      */
      virtual void getSolution();

      /**
      * Cleanup operation at the beginning of a sweep.
      */
      virtual void cleanup();

      /**
      * Has an association with the parent System been set?
      */
      bool hasSystem()
      {  return (systemPtr_ != 0); }

      /**
      * Return the parent system by reference.
      */
      System<D>& system()
      {  return *systemPtr_; }

      /// Whether to write real space concentration field files. 
      bool writeCRGrid_;

      /// Whether to write concentration field files in basis format. 
      bool writeCBasis_;

      /// Whether to write real space potential field files. 
      bool writeWRGrid_;

      // Protected members inherited from base classes
      using SweepTmpl< BasisFieldState<D> >::ns_;
      using SweepTmpl< BasisFieldState<D> >::baseFileName_;
      using SweepTmpl< BasisFieldState<D> >::initialize;
      using SweepTmpl< BasisFieldState<D> >::setCoefficients;
      using ParamComposite::readOptional;

   private:

      /// Trial state (produced by continuation in setGuess)
      BasisFieldState<D> trial_;

      /// Unit cell parameters for trial state 
      FSArray<double, 6> unitCellParameters_;

      /// Log file for summary output
      std::ofstream logFile_;

      /// Pointer to parent system.
      System<D>* systemPtr_;

      /// Output data to several files after convergence
      void outputSolution();

      /// Output brief summary of thermodynamic properties
      void outputSummary(std::ostream&);

   };

} // namespace Pspc
} // namespace Pscf
#endif
