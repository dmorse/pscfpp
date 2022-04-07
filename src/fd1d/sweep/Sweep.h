#ifndef FD1D_SWEEP_H
#define FD1D_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/sweep/SweepTmpl.h>             // base class
#include <fd1d/SystemAccess.h>                // base class
#include <fd1d/misc/HomogeneousComparison.h>  // member
#include <fd1d/misc/FieldIo.h>                // member
#include <util/containers/DArray.h>           // member

#include <util/global.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /**
   * Solve a sequence of problems along a line in parameter space.
   *
   * \ingroup Fd1d_Sweep_Module
   */
   class Sweep : public SweepTmpl< DArray<System::WField> >, 
                 public SystemAccess
   {

   public:

      typedef DArray<System::WField> State;

      typedef SweepTmpl<State>  Base;

      /**
      * Constructor.
      * 
      * \param system parent System object.
      */
      Sweep(System& system);

      /**
      * Destructor.
      */
      ~Sweep();

      /**
      * Read ns and baseFileName parameters.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Output information after obtaining a converged solution.
      *
      * \param stateFileName base name of output files
      */
      virtual void outputSolution(std::string const & stateFileName);

      /**
      * Output data to a running summary.
      *
      * \param outFile  output log file, open for writing
      */
      virtual void outputSummary(std::ostream& outFile);

   protected:

      /// Mode for comparison to homogeneous system (none -> -1)
      int homogeneousMode_;

      /**
      * Check allocation of w fields in one state, allocate if needed.
      *
      * \param state  one stored state of the system.
      */
      virtual void checkAllocation(State & state);

      /**
      * Setup operation at beginning sweep.
      */
      virtual void setup();

      /**
      * Set non-adjustable system parameters to new values.
      *
      * \param s new value of contour coordinate, in range [0,1]
      */
      virtual void setParameters(double s);

      /**
      * Create initial guess for new w fields by polynomial extrapolation.
      *
      * \param s new value of contour coordinate, in range [0,1]
      */
      virtual void extrapolate(double s);

      /**
      * Call the current iterator to solve one SCFT problem.
      *
      * Return 0 for sucessful solution, 1 on failure to converge.
      */
      virtual int solve(bool isContinuation);

      /**
      * Reset system to previous solution after iterature failure.
      *
      * The implementation of this function resets the system state to
      * correspond to that stored in state(0).
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
      * Close log file after end of sweep.
      */
      virtual void cleanup();

   private:

      /// Algorithm for comparing to a homogeneous system
      HomogeneousComparison comparison_;

      /// FieldIo object for writing output files
      FieldIo fieldIo_;

      /// Summary log file
      std::ofstream logFile_;

      /// Assign state rhs = lhs
      void assignFields(State& lhs, State const & rhs) const;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
