#ifndef RPC_AM_ITERATOR_BASIS_TPP
#define RPC_AM_ITERATOR_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/inter/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/containers/DArray.h>
#include <util/global.h>
#include <cmath>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
    : Iterator<D>(system)
   {
      isSymmetric_ = true;  
      ParamComposite::setClassName("AmIteratorBasis"); 
   }

   /*
   * Destructor.
   */
   template <int D>
   AmIteratorBasis<D>::~AmIteratorBasis()
   {}

   /*
   * Read parameters from file.
   */
   template <int D>
   void AmIteratorBasis<D>::readParameters(std::istream& in)
   {
      // Use base class methods to read parameters
      AmTmpl::readParameters(in);
      AmTmpl::readErrorType(in);

      UnitCell<D> const & unitCell = system().domain().unitCell();
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      int np = unitCell.nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);

      // Read optional isFlexible boolean (true by default)
      isFlexible_ = 1;  // Default
      readOptional(in, "isFlexible", isFlexible_);

      // Populate flexibleParams_ based on isFlexible_ (all 0s or all 1s),
      // then optionally overwrite with user input from param file
      if (isFlexible_) {
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(true); // Set all values to true
         }
         // Read optional flexibleParams_ array to overwrite current array
         readOptionalFSArray(in, "flexibleParams", flexibleParams_, np);
         if (nFlexibleParams() == 0) isFlexible_ = false;
      } else { // isFlexible_ = false
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(false); // Set all values to false
         }
      }

      // Read optional scaleStress value
      scaleStress_ = 10.0; // Default value
      readOptional(in, "scaleStress", scaleStress_);

      // Read option mixing parameters (lambda, useLambdaRamp, and r)
      AmTmpl::readMixingParameters(in);

      // Allocate local modified copy of Interaction class
      const int nMonomer = system().mixture().nMonomer();
      interaction_.setNMonomer(nMonomer);

   }

   /*
   * Output timing results to log file.
   */
   template<int D>
   void AmIteratorBasis<D>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "Iterator times contributions:\n";
      AmTmpl::outputTimers(out);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorBasis<D>::setup(bool isContinuation)
   {
      AmTmpl::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual functions to exchange data with parent system

   /*
   * Compute and return number of elements in a residual vector.
   */
   template <int D>
   int AmIteratorBasis<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();

      int nEle = nMonomer*nBasis;
      if (isFlexible()) {
         nEle += nFlexibleParams();
      }
      return nEle;
   }

   /*
   * Does the system have an initial field guess?
   */
   template <int D>
   bool AmIteratorBasis<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field vector (w fields and lattice parameters).
   */
   template <int D>
   void AmIteratorBasis<D>::getCurrent(DArray<double>& state)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();
      const int nEle = nElements();
      UTIL_CHECK(state.capacity() == nEle);

      // Copy w-fieldd into linear array
      int begin;
      for (int i = 0; i < nMonomer; i++) {
         DArray<double> const & field = system().w().basis(i);
         begin = i*nBasis;
         for (int k = 0; k < nBasis; k++) {
            state[begin + k] = field[k];
         }
      }

      // Add elements associated with unit cell parameters (if any)
      if (isFlexible()) {
         UTIL_CHECK(nFlexibleParams() > 0);
         UnitCell<D> const & unitCell = system().domain().unitCell();
         FSArray<double,6> const & parameters = unitCell.parameters();
         const int nParam = unitCell.nParameter();
         begin = nMonomer*nBasis;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               state[begin + counter] = scaleStress_*parameters[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
      }

   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void AmIteratorBasis<D>::evaluate()
   {  system().compute(isFlexible_); }

   /*
   * Compute the residual for the current system state.
   */
   template <int D>
   void AmIteratorBasis<D>::getResidual(DArray<double>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();
      const int n = nElements();
      UTIL_CHECK(resid.capacity() == n);

      // Local variables
      double chi, p;
      int i, j, k, begin;

      // Initialize residual vector to zero
      for (i = 0 ; i < n; ++i) {
         resid[i] = 0.0;
      }

      // Compute SCF residual vector elements
      for (i = 0; i < nMonomer; ++i) {
         begin = i*nBasis;
         for (j = 0; j < nMonomer; ++j) {
            chi = interaction_.chi(i,j);
            p = interaction_.p(i,j);
            DArray<double> const & c = system().c().basis(j);
            DArray<double> const & w = system().w().basis(j);
            if (system().h().hasData()) {
               DArray<double> const & h = system().h().basis(j);
               for (k = 0; k < nBasis; ++k) {
                  resid[begin + k] += chi*c[k] + p*(h[k] - w[k]);
               }
            } else {
               for (k = 0; k < nBasis; ++k) {
                  resid[begin + k] += chi*c[k] - p*w[k];
               }
            }
         }
      }

      // Add term proportional to sum of all concentrations
      double shift = -1.0 / interaction_.sumChiInverse();
      if (system().mask().hasData()) {
         DArray<double> const & mask = system().mask().basis();
         for (i = 0; i < nMonomer; ++i) {
            begin = i*nBasis;
            for (k = 0; k < nBasis; ++k) {
               resid[begin + k] += shift*mask[k];
            }
         }
      } else {
         for (i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] += shift;
         }
      }

      // If canonical ensemble, zero homogeneous residual components
      if (system().mixture().isCanonical()) {
         for (i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If flexible unit cell, then compute stress residuals
      if (isFlexible()) {

         // Combined -1 factor and stress scaling here. This is okay:
         // - residuals only show up as dot products (U, v, norm)
         //   or with their absolute value taken (max), so the
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale,
         //   and then undone right before passing to the unit cell.

         double coeff = -1.0 * scaleStress_;
         const int nParam = system().domain().unitCell().nParameter();
         begin = nMonomer*nBasis;
         int counter = 0;
         for (i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               resid[begin + counter] = coeff * stress(i);
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
      }

   }

   /*
   * Update the current system field vector.
   */
   template <int D>
   void AmIteratorBasis<D>::update(DArray<double>& newState)
   {
      // Constants
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();

      // Allocate wFields container
      DArray< DArray<double> > wFields;
      wFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         wFields[i].allocate(nBasis);
      }

      // Copy w fields from newState to wFields container
      int begin;
      for (int i = 0; i < nMonomer; i++) {
         begin = i*nBasis;
         for (int k = 0; k < nBasis; k++) {
            wFields[i][k] = newState[begin + k];
         }
      }

      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {

         // Set homogeneous components of all w fields to zero
         for (int i = 0; i < nMonomer; ++i) {
            wFields[i][0] = 0.0;
         }

         // Add average values arising from interactions
         double chi, wAve, cAve;
         for (int i = 0; i < nMonomer; ++i) {
            wAve = 0.0;
            for (int j = 0; j < nMonomer; ++j) {
               chi = interaction_.chi(i,j);
               cAve = system().c().basis(j)[0];
               wAve += chi * cAve;
            }
            wFields[i][0] = wAve;
         }

         // If external fields exist, add their spatial averages
         if (system().h().hasData()) {
            for (int i = 0; i < nMonomer; ++i) {
               wFields[i][0] += system().h().basis(i)[0];
            }
         }
      }

      // Set fields in system w container
      system().w().setBasis(wFields);

      // If flexible, update unit cell parameters
      if (isFlexible()) {

         // Initialize parameters array with current values
         FSArray<double, 6> parameters;
         parameters = system().domain().unitCell().parameters();

         // Reset any parameters that are flexible
         const double coeff = 1.0 / scaleStress_;
         const int nParam = system().domain().unitCell().nParameter();
         const int begin = nMonomer*nBasis;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = coeff * newState[begin + counter];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());

         // Set system unit cell parameters
         system().setUnitCell(parameters);
      }
   }

   /*
   * Output relevant system details to the iteration log.
   */
   template<int D>
   void AmIteratorBasis<D>::outputToLog()
   {
      if (isFlexible() && verbose() > 1) {
         double str;
         UnitCell<D> const & unitCell = system().domain().unitCell();
         const int nParam = unitCell.nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nBasis = system().domain().basis().nBasis();
         const int begin = nMonomer*nBasis;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               str = -1.0 * residual()[begin + counter] / scaleStress_;
               Log::file() 
                  << " Cell Param  " << i << " = "
                  << Dbl(unitCell.parameters()[i], 15)
                  << " , stress = " 
                  << Dbl(str, 15)
                  << "\n";
               counter++;
            }
         }
      }
   }

}
}
#endif
