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
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>
#include <pscf/inter/Interaction.h>
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
   * Read parameter file block.
   */
   template <int D>
   void AmIteratorBasis<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmTmpl::readParameters(in);
      AmTmpl::readErrorType(in);

      // Read optional isFlexible value
      isFlexible_ = 1; // default value
      readOptional(in, "isFlexible", isFlexible_);

      // Get and check the number of unit cell parameters
      int np = system().domain().unitCell().nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);
      UTIL_CHECK(system().domain().unitCell().lattice() 
                  != UnitCell<D>::Null);

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
      scaleStress_ = 10.0;  // default value
      readOptional(in, "scaleStress", scaleStress_);

      // Optionally read mixing parameters (lambda, useLambdaRamp, r)
      AmTmpl::readErrorType(in);

      // Optionally read mixing parameters (lambda, useLambdaRamp, r)
      AmTmpl::readMixingParameters(in);

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
      AmTmpl::outputTimers(out);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorBasis<D>::setup(bool isContinuation)
   {
      // Call parent setup method
      AmTmpl::setup(isContinuation);

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

      if (isFlexible_) {
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

      if (isFlexible_) {
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
         UTIL_CHECK(counter == nFlexibleParams());
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
      system().compute(isFlexible_);
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
      if (isFlexible_) {
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
               double str = stress(i);

               resid[nMonomer*nBasis + counter] = -1 * scaleStress_ * str;
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
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

      if (isFlexible_) {
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
         UTIL_CHECK(counter == nFlexibleParams());

         system().setUnitCell(parameters);
      }
   }

   /*
   * Output relevant system details to the iteration log file.
   */
   template<int D>
   void AmIteratorBasis<D>::outputToLog()
   {
      if (isFlexible_ && verbose() > 1) {
         const int nParam = system().domain().unitCell().nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nBasis = system().domain().basis().nBasis();
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               double str = residual()[nMonomer*nBasis + counter] /
                            (-1.0 * scaleStress_);
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

   #if 0
   // Private virtual functions for vector math

   /*
   * Set a vector equal to another (assign a = b).
   */
   template <int D>
   void AmIteratorBasis<D>::setEqual(DArray<double>& a,
                                     DArray<double> const & b)
   {  a = b; }

   /*
   * Compute the inner product of two real vectors.
   */
   template <int D>
   double AmIteratorBasis<D>::dotProduct(DArray<double> const & a,
                                         DArray<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = 0.0;
      for (int i=0; i < n; ++i) {
         // if either value is NaN, throw NanException
         if (std::isnan(a[i]) || std::isnan(b[i])) {
            throw NanException("AmIteratorBasis::dotProduct", __FILE__,
                               __LINE__, 0);
         }
         product += a[i] * b[i];
      }
      return product;
   }

   /*
   * Compute and return maximum element of residual vector.
   */
   template <int D>
   double AmIteratorBasis<D>::maxAbs(DArray<double> const & a)
   {
      const int n = a.capacity();
      double max = 0.0;
      double value;
      for (int i = 0; i < n; i++) {
         value = a[i];
         if (std::isnan(value)) { // if value is NaN, throw NanException
            throw NanException("AmIteratorBasis::dotProduct", __FILE__,
                               __LINE__, 0);
         }
         if (fabs(value) > max)
            max = fabs(value);
      }
      return max;
   }

   /*
   * Compute the vector difference a = b - c
   */
   template <int D>
   void AmIteratorBasis<D>::subVV(DArray<double>& a,
                                  DArray<double> const & b,
                                  DArray<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      UTIL_CHECK(n == c.capacity());
      for (int i = 0; i < n; i++) {
         a[i] = b[i] - c[i];
      }
   }

   /*
   * Composite a += b*c for vectors a and b, scalar c
   */
   template <int D>
   void AmIteratorBasis<D>::addEqVc(DArray<double>& a,
                                    DArray<double> const & b,
                                    double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      for (int i = 0; i < n; i++) {
         a[i] += c*b[i];
      }
   }
   #endif

}
}
#endif
