#ifndef RPG_AM_ITERATOR_BASIS_TPP
#define RPG_AM_ITERATOR_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include <rpg/System.h>

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

   // Constructor
   template <int D>
   AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
    : Iterator<D>(system),
      imposedFields_(system)
   {
      isSymmetric_ = true;  
      setClassName("AmIteratorBasis"); 
   }

   // Destructor
   template <int D>
   AmIteratorBasis<D>::~AmIteratorBasis()
   {}

   // Read parameter file block
   template <int D>
   void AmIteratorBasis<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters 
      AmIteratorTmpl<Iterator<D>,DArray<double> >::readParameters(in);
      AmIteratorTmpl<Iterator<D>,DArray<double> >::readErrorType(in);

      // Allocate local modified copy of Interaction class
      interaction_.setNMonomer(system().mixture().nMonomer());

      // Default parameter values
      isFlexible_ = 1;
      scaleStress_ = 10.0;

      int np = system().domain().unitCell().nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);
      UTIL_CHECK(system().domain().unitCell().lattice() != UnitCell<D>::Null);

      // Read in optional isFlexible value
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
      readOptional(in, "scaleStress", scaleStress_);

      // Read optional ImposedFieldsGenerator object
      readParamCompositeOptional(in, imposedFields_);
   }

   // Output timing results to log file.
   template<int D>
   void AmIteratorBasis<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Iterator times contributions:\n";
      AmIteratorTmpl<Iterator<D>, DArray<double> >::outputTimers(out);
   }

   // -- Protected virtual function -- //

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorBasis<D>::setup(bool isContinuation)
   {
      if (imposedFields_.isActive()) {
         imposedFields_.setup();
      }

      // Call parent setup method
      AmIteratorTmpl<Iterator<D>, DArray<double> >::setup(isContinuation);

      // Update chi matrix and related properties in member interaction_
      interaction_.update(system().interaction());
   }


   // -- Private virtual functions used to implement AM algorithm --  //

   // Set a vector equal to another (assign a = b)
   template <int D>
   void AmIteratorBasis<D>::setEqual(DArray<double>& a, 
                                     DArray<double> const & b)
   {  a = b; }

   // Compute the inner product of two real vectors.
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

   // Compute and return maximum element of residual vector.
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

   // Update the series of residual vectors.
   template <int D>
   void 
   AmIteratorBasis<D>::updateBasis(RingBuffer<DArray<double> > & basis,
                                   RingBuffer<DArray<double> > const& hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      DArray<double> newbasis;
      newbasis.allocate(n);

      // New basis vector is difference between two most recent states
      for (int i = 0; i < n; i++) {
         newbasis[i] = hists[0][i] - hists[1][i]; 
      }

      basis.append(newbasis);
   }

   // Compute trial field so as to minimize L2 norm of residual.
   template <int D>
   void
   AmIteratorBasis<D>::addHistories(DArray<double>& trial,
                                    RingBuffer<DArray<double> > const & basis,
                                    DArray<double> coeffs,
                                    int nHist)
   {
      int n = trial.capacity();
      for (int i = 0; i < nHist; i++) {
         for (int j = 0; j < n; j++) {
            // Not clear on the origin of the -1 factor
            trial[j] += coeffs[i] * -1 * basis[i][j];
         }
      }
   }

   // Add predicted error to the trial field.
   template <int D>
   void AmIteratorBasis<D>::addPredictedError(DArray<double>& fieldTrial,
                                              DArray<double> const & resTrial,
                                              double lambda)
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
   }

   // Does the system have an initial field guess?
   template <int D>
   bool AmIteratorBasis<D>::hasInitialGuess()
   {
      return system().w().hasData();
   }

   // Compute the number of elements in the residual vector.
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

   // Get the current w fields and lattice parameters
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

   // Perform the main system computation (solve the MDE)
   template <int D>
   void AmIteratorBasis<D>::evaluate()
   {
      // Solve MDEs for current omega field
      // (computes stress if isFlexible_ == true)
      system().compute(isFlexible_);
   }

   // Compute the residual for the current system state
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
      if (system().hasMask()) {
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
      if (system().hasExternalFields()) {
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
               double stress = system().mixture().stress(i);

               // Correct stress to account for effect of imposed fields
               if (imposedFields_.isActive()) {
                  stress = imposedFields_.correctedStress(i,stress);
               } 

               resid[nMonomer*nBasis + counter] = -1 * scaleStress_ * stress;
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
      }

   }

   // Update the current system field coordinates
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
         // If iterator has external fields, include them in homogeneous field
         if (system().hasExternalFields()) {
            for (int i = 0; i < nMonomer; ++i) {
               wField[i][0] += system().h().basis(i)[0];
            }
         }
      }
      system().setWBasis(wField);

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

      // Update imposed fields if needed
      if (imposedFields_.isActive()) {
         imposedFields_.update();
      }
   }

   // Output relevant system details to the iteration log file.
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
               double stress = residual()[nMonomer*nBasis + counter] /
                               (-1.0 * scaleStress_);
               Log::file() 
                      << " Cell Param  " << i << " = "
                      << Dbl(system().domain().unitCell().parameters()[i], 15)
                      << " , stress = " 
                      << Dbl(stress, 15)
                      << "\n";
               counter++;
            }
         }
      }
   }

   // Return specialized sweep parameter types to add to the Sweep object
   template<int D>
   GArray<ParameterType> AmIteratorBasis<D>::getParameterTypes()
   {
      GArray<ParameterType> arr;
      if (imposedFields_.isActive()) {
         arr = imposedFields_.getParameterTypes();
      } 
      return arr;
   }

   // Set the value of a specialized sweep parameter
   template<int D>
   void AmIteratorBasis<D>::setParameter(std::string name, DArray<int> ids, 
                                         double value, bool& success)
   {
      if (imposedFields_.isActive()) {
         imposedFields_.setParameter(name, ids, value, success);
      } else {
         success = false;
      }
   }

   // Get the value of a specialized sweep parameter
   template<int D>
   double AmIteratorBasis<D>::getParameter(std::string name, 
                                           DArray<int> ids, bool& success)
   const
   {
      if (imposedFields_.isActive()) {
         return imposedFields_.getParameter(name, ids, success);
      } else {
         success = false;
         return 0.0;
      }
   }

}
}
#endif
