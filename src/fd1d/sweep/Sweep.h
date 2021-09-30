#ifndef FD1D_SWEEP_H
#define FD1D_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class
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
   class Sweep : public ParamComposite, public SystemAccess
   {

   public:

      /**
      * Default Constructor.
      * 
      * Objects instantiated with this constructor must also call
      * the setSystem() function.
      */
      Sweep();

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
      * Set system if created with default constructor.
      */
      void setSystem(System& system);

      /**
      * Read ns and baseFileName parameters.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set system parameters.
      *
      * \param s path length coordinate, in range [0,1]
      */
      // virtual void setState(double s) = 0;

      /**
      * Output information after obtaining a converged solution.
      *
      * \param stateFileName base name of output files
      * \param s value of path length parameter s
      */
      virtual void outputSolution(std::string const & stateFileName, double s);

      /**
      * Output data to a running summary.
      *
      * \param outFile  output file, open for writing
      * \param i  integer index
      * \param s  value of path length parameter s
      */
      virtual void outputSummary(std::ostream& outFile, int i, double s);

      /**
      * Iterate to solution.
      */
      virtual void sweep();

   protected:

      /// Number of steps. 
      int ns_;

      /// Mode for comparison to homogeneous system (none -> -1)
      int homogeneousMode_;

      /// Base name for output files
      std::string baseFileName_;

      /**
      * Setup operation at beginning sweep.
      *
      * Must call initializeHistory.
      */
      virtual void setup();

      /**
      * Set non-adjustable system parameters to new values.
      *
      * \param s path length coordinate, in range [0,1]
      */
      virtual void setParameters(double s) = 0;

      /**
      * Create guess for adjustable variables by continuation.
      */
      virtual void setGuess(double s);

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

   private:

      /// Algorithm for comparing to a homogeneous system
      HomogeneousComparison comparison_;

      /// FieldIo object for writing output files
      FieldIo fieldIo_;

      /// Current solution 
      DArray<System::WField> wFields0_;

      /// Previous solution (for 1st order continuation)
      DArray<System::WField> wFields1_;

      void assignFields(DArray<System::WField>& lhs, 
                        DArray<System::WField> const & rhs) const;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
