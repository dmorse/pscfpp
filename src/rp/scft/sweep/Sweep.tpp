#ifndef RP_SWEEP_TPP
#define RP_SWEEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"                   // class template header

#include <rp/scft/iterator/Iterator.h>
#include <rpc/scft/ScftThermo.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/environment/Environment.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/interaction/Interaction.h>

#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <pscf/sweep/SweepTmpl.tpp>  // base class template implementation

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Maximum number of previous states = order of continuation + 1
   #define RP_HISTORY_CAPACITY 3

   /*
   * Default constructor (for unit testing).
   */
   template <int D, class T>
   Sweep<D,T>::Sweep()
    : SweepTmplT(RP_HISTORY_CAPACITY),
      writeCRGrid_(false),
      writeCBasis_(false),
      writeWRGrid_(false),
      systemPtr_(nullptr)
   {}

   /*
   * Constructor, creates association with parent system.
   */
   template <int D, class T>
   Sweep<D,T>::Sweep(System<D,T> & sys)
    : SweepTmplT(RP_HISTORY_CAPACITY),
      writeCRGrid_(false),
      writeCBasis_(false),
      writeWRGrid_(false),
      systemPtr_(&sys)
   {
      // Get specialized sweep parameters from Environment
      if (system().hasEnvironment()) {
         SweepTmplT::addParameterTypes(
                              system().environment().getParameterTypes());
      }
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   Sweep<D,T>::~Sweep()
   {}

   /*
   * Set association with a parent system (for unit testing).
   */
   template <int D, class T>
   void Sweep<D,T>::setSystem(System<D,T>& system)
   {  systemPtr_ = &system; }

   /*
   * Read parameters.
   */
   template <int D, class T>
   void Sweep<D,T>::readParameters(std::istream& in)
   {
      // Call the base class's readParameters function.
      SweepTmplT::readParameters(in);

      // Read optional flags indicating which field types to output
      ParamComposite::readOptional(in, "writeCRGrid", writeCRGrid_);
      ParamComposite::readOptional(in, "writeCBasis", writeCBasis_);
      ParamComposite::readOptional(in, "writeWRGrid", writeWRGrid_);
   }

   /*
   * Check allocation of one state object, allocate if necessary.
   */
   template <int D, class T>
   void Sweep<D,T>::checkAllocation(BasisFieldState<D,T>& state)
   {
      UTIL_CHECK(hasSystem());
      state.setSystem(system());
      state.allocate();
      state.unitCell() = system().domain().unitCell();
   }

   /*
   * Setup operations at the beginning of a sweep.
   */
   template <int D, class T>
   void Sweep<D,T>::setup()
   {
      SweepTmplT::initialize();
      checkAllocation(trial_);

      // Open log summary file
      std::string fileName = SweepTmplT::baseFileName_;
      fileName += "sweep.log";
      system().fileMaster().openOutputFile(fileName, logFile_);
      logFile_ << " step             ds     free_energy        pressure"
               << std::endl;
   };

   /*
   * Set non-adjustable system parameters to new values.
   *
   * \param s path length coordinate, in range [0,1]
   */
   template <int D, class T>
   void Sweep<D,T>::setParameters(double s)
   {
      // Empty default implementation allows Sweep<D> to be compiled.
      UTIL_THROW("Calling unimplemented function Sweep::setParameters");
   };

   /*
   * Create guess for adjustable variables by polynomial extrapolation.
   */
   template <int D, class T>
   void Sweep<D,T>::extrapolate(double sNew)
   {
      const int historySize = SweepTmplT::historySize();
      UTIL_CHECK(historySize > 0);

      // If historySize == 1, do nothing: Use previous system state
      // as trial for the new state.

      if (historySize > 1) {
         UTIL_CHECK(historySize <= SweepTmplT::historyCapacity());

         // Does the iterator allow a flexible unit cell ?
         bool isFlexible = system().iterator().isFlexible();

         // Compute coefficients of polynomial extrapolation to sNew
         SweepTmplT::setCoefficients(sNew);

         // Set extrapolated trial w fields
         double coeff;
         int nMonomer = system().mixture().nMonomer();
         int nBasis = system().domain().basis().nBasis();
         DArray<double>* newFieldPtr;
         DArray<double>* oldFieldPtr;
         int i, j, k;
         for (i=0; i < nMonomer; ++i) {
            newFieldPtr = &(trial_.field(i));

            // Previous state k = 0 (most recent)
            oldFieldPtr = &SweepTmplT::state(0).field(i);
            coeff = SweepTmplT::c(0);
            for (j=0; j < nBasis; ++j) {
               (*newFieldPtr)[j] = coeff*(*oldFieldPtr)[j];
            }

            // Previous states k >= 1 (older)
            for (k = 1; k < historySize; ++k) {
               oldFieldPtr = &SweepTmplT::state(k).field(i);
               coeff = SweepTmplT::c(k);
               for (j=0; j < nBasis; ++j) {
                  (*newFieldPtr)[j] += coeff*(*oldFieldPtr)[j];
               }
            }
         }

         // Make sure unitCellParameters_ is up to date with system
         // (if we are sweeping in a lattice parameter, then the system
         // parameters will be up-to-date but unitCellParameters_ wont be)
         FSArray<double, 6> oldParameters = unitCellParameters_;
         unitCellParameters_ = system().domain().unitCell().parameters();

         // If isFlexible, extrapolate the flexible unit cell parameters
         if (isFlexible) {

            double coeff;
            double parameter;
            const FSArray<bool,6> flexParams 
                            = system().iterator().flexibleParams();
            const int nParameter
                            = system().domain().unitCell().nParameter();
            UTIL_CHECK(flexParams.size() == nParameter);

            // Add contributions from all previous states
            for (k = 0; k < historySize; ++k) {
               coeff = SweepTmplT::c(k);
               for (int i = 0; i < nParameter; ++i) {
                  if (flexParams[i]) {
                     if (k == 0) {
                        unitCellParameters_[i] = 0.0;
                     }
                     parameter 
                          = SweepTmplT::state(k).unitCell().parameter(i);
                     unitCellParameters_[i] += coeff*parameter;
                  }
               }
            }

         }

         // Reset trial_.unitCell() object
         trial_.unitCell().setParameters(unitCellParameters_);

         // Check if any unit cell parameters have changed
         bool newCellParams(true);
         for (int i = 0; i < oldParameters.size(); i++) {
            if (fabs(oldParameters[i] - unitCellParameters_[i]) < 1e-10) {
               newCellParams = false;
               break;
            }
         }

         // Transfer data from trial_ state to parent system
         trial_.setSystemState(newCellParams);
      }

   };

   /*
   * Call current iterator to solve SCFT problem.
   *
   * Return 0 for sucessful solution, 1 on failure to converge.
   */
   template <int D, class T>
   int Sweep<D,T>::solve(bool isContinuation)
   {  return system().iterate(isContinuation); };

   /*
   * Reset system to previous solution after iterature failure.
   *
   * The implementation of this function should reset the system state
   * to correspond to that stored in state(0).
   */
   template <int D, class T>
   void Sweep<D,T>::reset()
   {
      bool isFlexible = system().iterator().isFlexible();
      SweepTmplT::state(0).setSystemState(isFlexible);
   }

   /*
   * Update state(0) and output data after successful convergence
   *
   * The implementation of this function should copy the current
   * system state into state(0) and output any desired information
   * about the current converged solution.
   */
   template <int D, class T>
   void Sweep<D,T>::getSolution()
   {
      SweepTmplT::state(0).setSystem(system());
      SweepTmplT::state(0).getSystemState();

      // Output converged solution to several files
      outputSolution();

      // Output summary to log file
      outputSummary(logFile_);

   };

   template <int D, class T>
   void Sweep<D,T>::outputSolution()
   {
      std::ofstream out;
      std::string outFileName;
      std::string indexString = toString(SweepTmplT::nAccept() - 1);

      // Open parameter file, with thermodynamic properties at end
      outFileName = SweepTmplT::baseFileName_;
      outFileName += indexString;
      outFileName += ".stt";
      system().fileMaster().openOutputFile(outFileName, out);

      // Write data file, with thermodynamic properties at end
      system().writeParamNoSweep(out);
      out << std::endl;
      system().scft().write(out);
      out.close();

      // Write w fields
      UTIL_CHECK(system().w().hasData());
      outFileName = SweepTmplT::baseFileName_;
      outFileName += indexString;
      outFileName += "_w";
      if (system().w().isSymmetric()) {
         outFileName += ".bf";
         system().w().writeBasis(outFileName);
      } else {
         outFileName += ".rf";
         system().w().writeRGrid(outFileName);
      }

      // Optionally write c rgrid files
      if (writeCRGrid_) {
         UTIL_CHECK(system().c().hasData());
         outFileName = SweepTmplT::baseFileName_;
         outFileName += indexString;
         outFileName += "_c";
         outFileName += ".rf";
         system().c().writeRGrid(outFileName);
      }

      // Optionally write c basis files
      if (writeCBasis_ && system().c().isSymmetric()) {
         UTIL_CHECK(system().c().hasData());
         outFileName = SweepTmplT::baseFileName_;
         outFileName += indexString;
         outFileName += "_c";
         outFileName += ".bf";
         system().c().writeBasis(outFileName);
      }

      // Optionally write w rgrid files
      if (writeWRGrid_ && system().w().isSymmetric()) {
         outFileName = SweepTmplT::baseFileName_;
         outFileName += indexString;
         outFileName += "_w";
         outFileName += ".rf";
         system().w().writeRGrid(outFileName);
      }

   }

   template <int D, class T>
   void Sweep<D,T>::outputSummary(std::ostream& out)
   {
      int i = SweepTmplT::nAccept() - 1;
      double sNew = SweepTmplT::s(0);
      if (!system().scft().hasData()) system().scft().compute();
      out << Int(i,5) << Dbl(sNew)
          << Dbl(system().scft().fHelmholtz(),16)
          << Dbl(system().scft().pressure(),16);
      out << std::endl;
   }

   template <int D, class T>
   void Sweep<D,T>::cleanup()
   {  logFile_.close(); }

} // namespace Rp
} // namespace Pscf
#endif
