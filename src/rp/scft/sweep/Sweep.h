#ifndef RP_SWEEP_H
#define RP_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/sweep/SweepTmpl.h>          // base class template
#include <util/containers/FSArray.h>       // member
#include <util/global.h>                   // inline functions

#include <fstream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Solve a sequence of SCFT problems along a line in parameter space.
   *
   * Instantiations of this class template are used as base classes for 
   * two closely analogous class templates, also named Sweep, that are
   * defined in Rpc and Rpg namespaces and used in the pscf_rpc and 
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref scft_sweep_page "Manual page"
   * \ingroup Rp_Scft_Sweep_Module
   */
   template <int D, class T>
   class Sweep : public SweepTmpl< typename T::BasisFieldState >
   {

   public:

      /**
      * Set association with parent system.
      *
      * Call for objects created with default constructor.
      *
      * \param system  parent system
      */
      void setSystem(typename T::System& system);

      /**
      * Read parameters from param file.
      *
      * \param in  parameter file input stream
      */
      virtual void readParameters(std::istream& in);

      // Public members inherited from base class template SweepTmpl
      using SweepTmpl< typename T::BasisFieldState >::addParameterTypes;
      using SweepTmpl< typename T::BasisFieldState >::addParameterType;
      using SweepTmpl< typename T::BasisFieldState >::historyCapacity;
      using SweepTmpl< typename T::BasisFieldState >::historySize;
      using SweepTmpl< typename T::BasisFieldState >::nAccept;
      using SweepTmpl< typename T::BasisFieldState >::state;
      using SweepTmpl< typename T::BasisFieldState >::s;
      using SweepTmpl< typename T::BasisFieldState >::c;

   protected:

      /**
      * Default constructor.
      */
      Sweep();

      /**
      * Constructor, creates assocation with parent system.
      *
      * \param system  parent system
      */
      Sweep(typename T::System& system);

      /**
      * Destructor.
      */
      ~Sweep();

      /**
      * Check allocation of fields in one state, allocate if necessary.
      *
      * \param state  stored state of the system
      */
      virtual void checkAllocation(typename T::BasisFieldState& state);

      /**
      * Setup operation at the beginning of a sweep.
      */
      virtual void setup();

      /**
      * Set system parameters to new values.
      *
      * \param sNew contour variable value for new trial solution.
      */
      virtual void setParameters(double sNew) = 0;

      /**
      * Create a guess for adjustable variables by continuation.
      *
      * The "adjustable variables" in a standard SCFT problem are w field
      * values and adjustable unit cell parameters, i.e., variables that
      * are adjusted by the iterator.
      *
      * \param sNew  contour variable value for new trial solution.
      */
      virtual void extrapolate(double sNew);

      /**
      * Call current iterator to solve SCFT problem.
      *
      * \param isContinuation  true iff is continuation in a sweep
      * \return 0 for sucessful solution, 1 on failure to converge
      */
      virtual int solve(bool isContinuation);

      /**
      * Reset system to previous solution after iterature failure.
      *
      * The implementation of this function should reset the system state
      * to correspond to that stored in state(0), i.e., the previous
      * converged solution.
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
      * Does an association with the parent system exist?
      */
      bool hasSystem()
      {  return (bool)(systemPtr_); }

      /**
      * Return the parent system by reference.
      */
      typename T::System& system()
      {
         UTIL_CHECK(systemPtr_);
         return *systemPtr_;
      }

      // Protected variables writeCGrid_, writeCBasis_, and writeWGrid_
      // control which converged fields will be written to files after
      // solution of each SCFT problem within a sweep.

      /**
      * Should concentration fields be written to file in r-grid format?
      */
      bool writeCRGrid_;

      /**
      * Should concentration fields be written to file in basis format?
      */
      bool writeCBasis_;

      /**
      * Should converged w fields be written to file in r-grid format?
      */
      bool writeWRGrid_;

      // Protected inherited member variables
      using SweepTmplT = SweepTmpl< typename T::BasisFieldState >;
      using SweepTmplT::ns_;
      using SweepTmplT::baseFileName_;
      using SweepTmplT::initialize;
      using SweepTmplT::setCoefficients;
      using ParamComposite::readOptional;

   private:

      /// Trial state (produced by continuation in setGuess)
      typename T::BasisFieldState trial_;

      /// Unit cell parameters for trial state
      FSArray<double, 6> unitCellParameters_;

      /// Log file for summary output
      std::ofstream logFile_;

      /// Pointer to parent system
      typename T::System* systemPtr_;

      /// Output data to several files after convergence.
      void outputSolution();

      /// Output summary of thermodynamic properties.
      void outputSummary(std::ostream&);

   };

} // namespace Rp
} // namespace Pscf
#endif
