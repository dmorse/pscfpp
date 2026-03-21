#ifndef RPG_AM_ITERATOR_BASIS_TPP
#define RPG_AM_ITERATOR_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/global.h>
#include <cmath>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
    : AmIteratorTmplT(),
      interaction_(),
      scaleStress_(1.0)
   {
      Iterator<D>::setSystem(system);
      Iterator<D>::isSymmetric_ = true;
      ParamComposite::setClassName("AmIteratorBasis");
   }

   /*
   * Destructor.
   */
   template <int D>
   AmIteratorBasis<D>::~AmIteratorBasis()
   {}

   /*
   * Read parameter file block.
   */
   template <int D>
   void AmIteratorBasis<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmplT::readParameters(in);
      AmIteratorTmplT::readErrorType(in);

      // Read optional isFlexible value
      Iterator<D>::isFlexible_ = 1; // default value
      readOptional(in, "isFlexible", Iterator<D>::isFlexible_);

      // Get and check the number of unit cell parameters
      int np = system().domain().unitCell().nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);
      UTIL_CHECK(system().domain().unitCell().lattice() 
                  != UnitCell<D>::Null);

      // Populate flexibleParams_ based on isFlexible_ (all 0s or all 1s),
      // then optionally overwrite with user input from param file
      if (Iterator<D>::isFlexible_) {
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(true); // Set all values to true
         }
         // Read optional flexibleParams_ array to overwrite current array
         readOptionalFSArray(in, "flexibleParams", flexibleParams_, np);
         if (Iterator<D>::nFlexibleParams() == 0) {
            Iterator<D>::isFlexible_ = false;
         }
      } else { // isFlexible_ = false
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(false); // Set all values to false
         }
      }

      // Read optional scaleStress value
      scaleStress_ = 10.0;  // default value
      readOptional(in, "scaleStress", scaleStress_);

      // Optionally read mixing parameters (lambda, useLambdaRamp, r)
      AmIteratorTmplT::readErrorType(in);

      // Optionally read mixing parameters (lambda, useLambdaRamp, r)
      AmIteratorTmplT::readMixingParameters(in);

      // Allocate local modified copy of Interaction class
      interaction_.setNMonomer(system().mixture().nMonomer());
   }

   /*
   * Output timing results to log file.
   */
   template<int D>
   void AmIteratorBasis<D>::outputTimers(std::ostream& out) const
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Iterator times contributions:\n";
      AmIteratorTmplT::outputTimers(out);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorBasis<D>::setup(bool isContinuation)
   {
      // Call parent setup method
      AmIteratorTmplT::setup(isContinuation);

      // Update chi matrix and related properties in member interaction_
      interaction_.update(system().interaction());
   }

   // Private virtual functions that interact with parent system

   /*
   * Compute the number of elements in the residual vector.
   */
   template <int D>
   int AmIteratorBasis<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();

      int nEle = nMonomer*nBasis;

      if (Iterator<D>::isFlexible_) {
         nEle += Iterator<D>::nFlexibleParams();
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
   * Get the current w fields and lattice parameters.
   */
   template <int D>
   void AmIteratorBasis<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();
      const DArray< DArray<double> > & currSys = system().w().basis();

      for (int i = 0; i < nMonomer; i++) {
         for (int k = 0; k < nBasis; k++) {
            curr[i*nBasis+k] = currSys[i][k];
         }
      }

      if (Iterator<D>::isFlexible_) {
         const int begin = nMonomer*nBasis;
         const int nParam = system().domain().unitCell().nParameter();
         FSArray<double,6> const & parameters
                               = system().domain().unitCell().parameters();
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               curr[begin + counter] = scaleStress_ * parameters[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == Iterator<D>::nFlexibleParams());
      }

   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void AmIteratorBasis<D>::evaluate()
   {
      // Solve MDEs for current omega field
      // (computes stress if isFlexible_ == true)
      system().compute(Iterator<D>::isFlexible_);
   }

   /*
   * Compute the residual for the current system state.
   */
   template <int D>
   void AmIteratorBasis<D>::getResidual(DArray<double>& resid)
   {
      UTIL_CHECK(system().domain().basis().isInitialized());
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();
      const int n = nElements();

      // Initialize residuals
      for (int i = 0 ; i < n; ++i) {
         resid[i] = 0.0;
      }

      // Compute SCF residual vector elements
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            double chi = interaction_.chi(i,j);
            double p = interaction_.p(i,j);
            DArray<double> const & c = system().c().basis(j);
            DArray<double> const & w = system().w().basis(j);
            for (int k = 0; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] += chi*c[k] - p*w[k];
            }
         }
      }

      // If iterator has mask, account for it in residual values
      if (system().mask().hasData()) {
         DArray<double> const & mask = system().mask().basis();
         double sumChiInv = interaction_.sumChiInverse();
         for (int i = 0; i < nMonomer; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] -= mask[k] / sumChiInv;
            }
         }
      }

      // If iterator has external fields, account for them in the values
      // of the residuals
      if (system().h().hasData()) {
         for (int i = 0; i < nMonomer; ++i) {
            for (int j = 0; j < nMonomer; ++j) {
               double p = interaction_.p(i,j);
               DArray<double> const & h = system().h().basis(j);
               for (int k = 0; k < nBasis; ++k) {
                  int idx = i*nBasis + k;
                  resid[idx] += p * h[k];
               }
            }
         }
      }

      // If not canonical, account for incompressibility
      if (!system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= 1.0 / interaction_.sumChiInverse();
         }
      } else {
         // Explicitly set homogeneous residual components
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (Iterator<D>::isFlexible_) {
         const int nParam = system().domain().unitCell().nParameter();

         //  Note:
         //  Combined -1 factor and stress scaling here.  This is okay:
         //  - Residuals only show up as dot products (U, v, norm)
         //    or with their absolute value taken (max), so the
         //    sign on a given residual vector element is not relevant
         //    as long as it is consistent across all vectors
         //  - The scaling is applied here and to the unit cell param
         //    storage, so that updating is done on the same scale,
         //    and then undone right before passing to the unit cell.

         int counter = 0;
         for (int i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               double str = Iterator<D>::stress(i);
               resid[nMonomer*nBasis + counter] = -1 * scaleStress_ * str;
               counter++;
            }
         }
         UTIL_CHECK(counter == Iterator<D>::nFlexibleParams());
      }

   }

   /*
   * Update the current system field coordinates.
   */
   template <int D>
   void AmIteratorBasis<D>::update(DArray<double>& newGuess)
   {
      UTIL_CHECK(system().domain().basis().isInitialized());
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();

      DArray< DArray<double> > wField;
      wField.allocate(nMonomer);

      // Restructure in format of monomers, basis functions
      for (int i = 0; i < nMonomer; i++) {
         wField[i].allocate(nBasis);
         for (int k = 0; k < nBasis; k++) {
            wField[i][k] = newGuess[i*nBasis + k];
         }
      }

      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {
         double chi;
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               chi = interaction_.chi(i,j);
               wField[i][0] += chi * system().c().basis(j)[0];
            }
         }
         // If system has external fields, include them in homogeneous part
         if (system().h().hasData()) {
            for (int i = 0; i < nMonomer; ++i) {
               wField[i][0] += system().h().basis(i)[0];
            }
         }
      }
      system().w().setBasis(wField);

      if (Iterator<D>::isFlexible_) {
         const int nParam = system().domain().unitCell().nParameter();
         const int begin = nMonomer*nBasis;

         FSArray<double,6> parameters;
         parameters = system().domain().unitCell().parameters();

         const double coeff = 1.0 / scaleStress_;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = coeff * newGuess[begin + counter];
               counter++;
            }
         }
         UTIL_CHECK(counter == Iterator<D>::nFlexibleParams());

         system().setUnitCell(parameters);
      }
   }

   /*
   * Output relevant system details to the iteration log file.
   */
   template<int D>
   void AmIteratorBasis<D>::outputToLog()
   {
      if (Iterator<D>::isFlexible_ && AmIteratorTmplT::verbose() > 1) {
         const int nParam = system().domain().unitCell().nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nBasis = system().domain().basis().nBasis();
         double res, str;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               res = AmIteratorTmplT::residual()[nMonomer*nBasis + counter];
               str = -1.0 * res / scaleStress_;
               Log::file()
                   << " Cell Param  " << i << " = "
                   << Dbl(system().domain().unitCell().parameters()[i], 15)
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
