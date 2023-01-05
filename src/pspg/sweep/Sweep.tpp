#ifndef PSPG_SWEEP_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include <pscf/sweep/SweepTmpl.tpp>
#include <pspg/System.h>
#include <pspg/iterator/Iterator.h>
#include <pscf/inter/Interaction.h>

#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Maximum number of previous states = order of continuation + 1
   #define PSPG_HISTORY_CAPACITY 3

   /*
   * Default constructor.
   */
   template <int D>
   Sweep<D>::Sweep() 
    : SweepTmpl< BasisFieldState<D> >(PSPG_HISTORY_CAPACITY),
      systemPtr_(0)
   {}

   /*
   * Constructor, creates association with parent system.
   */
   template <int D>
   Sweep<D>::Sweep(System<D> & system) 
    : SweepTmpl< BasisFieldState<D> >(PSPG_HISTORY_CAPACITY),
      systemPtr_(&system)
   {}

   /*
   * Destructor.
   */
   template <int D>
   Sweep<D>::~Sweep() 
   {}

   /*
   * Set association with a parent system.
   */
   template <int D>
   void Sweep<D>::setSystem(System<D>& system) 
   {  systemPtr_ = &system; }

   /*
   * Check allocation of one state object, allocate if necessary.
   */
   template <int D>
   void Sweep<D>::checkAllocation(BasisFieldState<D>& state) 
   { 
      UTIL_CHECK(hasSystem()); 
      state.setSystem(system());
      state.allocate(); 
      state.unitCell() = system().unitCell();
   }

   /*
   * Setup operations at the beginning of a sweep.
   */
   template <int D>
   void Sweep<D>::setup() 
   {
      initialize();
      checkAllocation(trial_);

      // Open log summary file
      std::string fileName = baseFileName_;
      fileName += "sweep.log";
      system().fileMaster().openOutputFile(fileName, logFile_);
   };

   /*
   * Set non-adjustable system parameters to new values.
   *
   * \param s path length coordinate, in range [0,1]
   */
   template <int D>
   void Sweep<D>::setParameters(double s) 
   {
      // Empty default implementation allows Sweep<D> to be compiled.
      UTIL_THROW("Calling unimplemented function Sweep::setParameters");
   };

   /*
   * Create guess for adjustable variables by polynomial extrapolation.
   */
   template <int D>
   void Sweep<D>::extrapolate(double sNew) 
   {
      UTIL_CHECK(historySize() > 0);

      // If historySize() == 1, do nothing: Use previous system state
      // as trial for the new state.
      
      if (historySize() > 1) {
         UTIL_CHECK(historySize() <= historyCapacity());

         // Does the iterator allow a flexible unit cell ?
         bool isFlexible = system().iterator().isFlexible();

         // Compute coefficients of polynomial extrapolation to sNew
         setCoefficients(sNew);

         // Set extrapolated trial w fields 
         double coeff;
         int nMonomer = system().mixture().nMonomer();
         int nBasis = system().basis().nBasis();
         DArray<double>* newFieldPtr;
         DArray<double>* oldFieldPtr; 
         int i, j, k;
         for (i=0; i < nMonomer; ++i) {
            newFieldPtr = &(trial_.field(i));

            // Previous state k = 0 (most recent)
            oldFieldPtr = &state(0).field(i);
            coeff = c(0);
            for (j=0; j < nBasis; ++j) {
               (*newFieldPtr)[j] = coeff*(*oldFieldPtr)[j];
            }

            // Previous states k >= 1 (older)
            for (k = 1; k < historySize(); ++k) {
               oldFieldPtr = &state(k).field(i);
               coeff = c(k);
               for (j=0; j < nBasis; ++j) {
                  (*newFieldPtr)[j] += coeff*(*oldFieldPtr)[j];
               }
            }
         }

         // Make sure unitCellParameters_ is up to date with system
         // (if we are sweeping in a lattice parameter, then the system
         // parameters will be up-to-date but unitCellParameters_ wont be)
         FSArray<double, 6> oldParameters = unitCellParameters_;
         unitCellParameters_ = system().unitCell().parameters();

         // If isFlexible, then extrapolate the flexible unit cell parameters
         if (isFlexible) {

            double coeff;
            double parameter;
           
            #if 0 
            const FSArray<int,6> indices = system().iterator().flexibleParams();
            const int nParameter = indices.size();

            // Add contributions from all previous states
            for (k = 0; k < historySize(); ++k) {
               coeff = c(k);
               for (int i = 0; i < nParameter; ++i) {
                  if (k == 0) {
                     unitCellParameters_[indices[i]] = 0;
                  }
                  parameter = state(k).unitCell().parameter(indices[i]);
                  unitCellParameters_[indices[i]] += coeff*parameter;
               }
            }
            #endif

            // Add contributions from all previous states
            const int nParameter = unitCellParameters_.size();
            for (k = 0; k < historySize(); ++k) {
               coeff = c(k);
               for (int i = 0; i < nParameter; ++i) {
                  if (k == 0) {
                     unitCellParameters_[i] = 0.0;
                  }
                  parameter = state(k).unitCell().parameter(i);
                  unitCellParameters_[i] += coeff*parameter;
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
   template <int D>
   int Sweep<D>::solve(bool isContinuation)
   {  return system().iterate(isContinuation); };

   /*
   * Reset system to previous solution after iterature failure.
   *
   * The implementation of this function should reset the system state
   * to correspond to that stored in state(0).
   */
   template <int D>
   void Sweep<D>::reset()
   {  
      bool isFlexible = system().iterator().isFlexible();
      state(0).setSystemState(isFlexible); 
   }

   /*
   * Update state(0) and output data after successful convergence
   *
   * The implementation of this function should copy the current 
   * system state into state(0) and output any desired information
   * about the current converged solution.
   */
   template <int D>
   void Sweep<D>::getSolution() 
   { 
      state(0).setSystem(system());
      state(0).getSystemState(); 

      // Output converged solution to several files
      outputSolution();

      // Output summary to log file
      outputSummary(logFile_);

   };

   template <int D>
   void Sweep<D>::outputSolution()
   {
      std::ofstream out;
      std::string outFileName;
      std::string indexString = toString(nAccept() - 1);

      // Open parameter file, with thermodynamic properties at end
      outFileName = baseFileName_;
      outFileName += indexString;
      outFileName += ".dat";
      system().fileMaster().openOutputFile(outFileName, out);

      // Write data file, with thermodynamic properties at end
      out << "System{" << std::endl;
      system().mixture().writeParam(out);
      system().interaction().writeParam(out);
      out << "}" << std::endl;
      out << std::endl;
      out << "unitCell       " << system().unitCell();
      system().writeThermo(out);
      out.close();

      // Write w fields
      outFileName = baseFileName_;
      outFileName += indexString;
      outFileName += "_w";
      outFileName += ".bf";
      system().writeWBasis(outFileName);

      // Optionally write c rgrid files
      if (writeCRGrid_) {
        outFileName = baseFileName_;
        outFileName += indexString;
        outFileName += "_c";
        outFileName += ".rf";
        system().writeCRGrid(outFileName);
      }

      // Optionally write c basis files
      if (writeCBasis_) {
         outFileName = baseFileName_;
         outFileName += indexString;
         outFileName += "_c";
         outFileName += ".bf";
         system().writeCBasis(outFileName);
      }

      // Optionally write w rgrid files
      if (writeWRGrid_) {
        outFileName = baseFileName_;
        outFileName += indexString;
        outFileName += "_w";
        outFileName += ".rf";
        system().writeWRGrid(outFileName);
      }
   }

   template <int D>
   void Sweep<D>::outputSummary(std::ostream& out)
   {
      int i = nAccept() - 1;
      double sNew = s(0);
      out << Int(i,5) << Dbl(sNew)
          << Dbl(system().fHelmholtz(),16)
          << Dbl(system().pressure(),16);
      out << std::endl;
   }

   template <int D>
   void Sweep<D>::cleanup() 
   {  logFile_.close(); }

} // namespace Pspg
} // namespace Pscf
#endif
