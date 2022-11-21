#ifndef PSPC_AM_ITERATOR_TPP
#define PSPC_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <pspc/System.h>
#include <pscf/inter/Interaction.h>
#include <util/global.h>

namespace Pscf{
namespace Pspc{

   using namespace Util;

   // Constructor
   template <int D>
   AmIterator<D>::AmIterator(System<D>& system)
    : Iterator<D>(system)
   {  setClassName("AmIterator"); }

   // Destructor
   template <int D>
   AmIterator<D>::~AmIterator()
   {  }

   // Read parameters from file
   template <int D>
   void AmIterator<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Iterator<D>, DArray<double> >::readParameters(in);
      AmIteratorTmpl<Iterator<D>, DArray<double> >::readErrorType(in);

      // Default parameter values
      isFlexible_ = 0; 
      scaleStress_ = 10.0;

      int np = system().unitCell().nParameter();

      // Read optional isFlexible boolean
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
   }

   // Private virtual functions used to implement AM algorithm

   template <int D>
   void AmIterator<D>::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

   // Compute and return inner product of two vectors.
   template <int D>
   double AmIterator<D>::dotProduct(DArray<double> const & a, 
                                    DArray<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = 0.0;
      for (int i = 0; i < n; i++) {
         product += a[i] * b[i];
      }
      return product;
   }


   // Compute and return maximum element of a vector.
   template <int D>
   double AmIterator<D>::maxAbs(DArray<double> const & hist)
   {
      const int n = hist.capacity();
      double max = 0.0;
      for (int i = 0; i < n; i++) {
         if (fabs(hist[i]) > max)
            max = fabs(hist[i]);
      }
      return max;
   }

   #if 0
   // Compute and return L2 norm of residual vector
   template <int D>
   double AmIterator<D>::findNorm(DArray<double> const & hist)
   {
      const int n = hist.capacity();
      double normResSq = 0.0;

      for (int i = 0; i < n; i++) {
         normResSq += hist[i] * hist[i];
      }

      return sqrt(normResSq);
   }

   // Compute one element of U matrix of by computing a dot product
   template <int D>
   double
   AmIterator<D>::computeUDotProd(RingBuffer<DArray<double> > const & resBasis,
                                  int m, int n)
   {
      const int length = resBasis[0].capacity();

      double dotprod = 0.0;
      for(int i = 0; i < length; i++) {
         dotprod += resBasis[m][i] * resBasis[n][i];
      }

      return dotprod;
   }

   // Compute one element of V vector by computing a dot product
   template <int D>
   double
   AmIterator<D>::computeVDotProd(DArray<double> const & resCurrent,
                                  RingBuffer<DArray<double> > const & resBasis,
                                  int m)
   {
      const int length = resBasis[0].capacity();

      double dotprod = 0.0;
      for(int i = 0; i < length; i++) {
         dotprod += resCurrent[i] * resBasis[m][i];
      }

      return dotprod;
   }

   // Update entire U matrix
   template <int D>
   void AmIterator<D>::updateU(DMatrix<double> & U,
                               RingBuffer<DArray<double> > const & resBasis,
                               int nHist)
   {
      // Update matrix U by shifting elements diagonally
      int maxHist = U.capacity1();
      for (int m = maxHist-1; m > 0; --m) {
         for (int n = maxHist-1; n > 0; --n) {
            U(m,n) = U(m-1,n-1);
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nHist; ++m) {
         double dotprod = computeUDotProd(resBasis,0,m);
         U(m,0) = dotprod;
         U(0,m) = dotprod;
      }
   }

   template <int D>
   void AmIterator<D>::updateV(DArray<double> & v,
                               DArray<double> const & resCurrent,
                               RingBuffer<DArray<double> > const & resBasis,
                               int nHist)
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }
   #endif

   // Update basis
   template <int D>
   void AmIterator<D>::updateBasis(RingBuffer<DArray<double> > & basis,
                                   RingBuffer<DArray<double> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      DArray<double> newbasis;
      newbasis.allocate(n);

      for (int i = 0; i < n; i++) {
         // sequential histories basis vectors
         newbasis[i] = hists[0][i] - hists[1][i]; 
      }

      basis.append(newbasis);
   }

   template <int D>
   void
   AmIterator<D>::addHistories(DArray<double>& trial,
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

   template <int D>
   void AmIterator<D>::addPredictedError(DArray<double>& fieldTrial,
                                         DArray<double> const & resTrial,
                                         double lambda)
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
   }

   // Private virtual functions to exchange data with parent system

   // Does the system have an initial field guess?
   template <int D>
   bool AmIterator<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return number of elements in a residual vector
   template <int D>
   int AmIterator<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      int nEle = nMonomer*nBasis;

      if (isFlexible()) {
         nEle += nFlexibleParams();
      }

      return nEle;
   }

   // Get the current field from the system
   template <int D>
   void AmIterator<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      const DArray< DArray<double> > * currSys = &system().w().basis();

      for (int i = 0; i < nMonomer; i++) {
         for (int k = 0; k < nBasis; k++)
         {
            curr[i*nBasis+k] = (*currSys)[i][k];
         }
      }

      const int nParam = system().unitCell().nParameter();
      const FSArray<double,6> currParam = system().unitCell().parameters();

      int counter = 0;
      for (int i = 0; i < nParam; i++) {
         if (flexibleParams_[i]) {
            curr[nMonomer*nBasis + counter] = scaleStress_*currParam[i];
            counter++;
         }
      }
      UTIL_CHECK(counter == nFlexibleParams());

      return;
   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void AmIterator<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute();

      // Compute stress if required
      if (isFlexible()) {
         system().mixture().computeStress();
      }
   }

   // Compute the residual for the current system state
   template <int D>
   void AmIterator<D>::getResidual(DArray<double>& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      // Initialize residuals
      for (int i = 0 ; i < n; ++i) {
         resid[i] = 0.0;
      }

      // Compute SCF residual vector elements
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = 0; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] +=
                  system().interaction().chi(i,j)*system().c().basis(j)[k] -
                  system().interaction().idemp(i,j)*system().w().basis(j)[k];
            }
         }
      }

      // If iterator has mask, account for it in residual values
      if (system().hasMask()) {
         for (int i = 0; i < nMonomer; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] -= system().mask().basis()[k] / 
                             system().interaction().sum_inv();
            }
         }
      }

      // If iterator has external fields, account for them in the values 
      // of the residuals
      if (system().hasExternalFields()) {
         for (int i = 0; i < nMonomer; ++i) {
            for (int j = 0; j < nMonomer; ++j) {
               for (int k = 0; k < nBasis; ++k) {
                  int idx = i*nBasis + k;
                  resid[idx] += system().interaction().idemp(i,j) * 
                                system().h().basis(j)[k];
               }
            }
         }
      }

      // If not canonical, account for incompressibility
      if (!system().mixture().isCanonical()) {
         // Fraction of unit cell occupied by polymers (1 if no mask present)
         double phiTot = system().mask().phiTot();

         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= phiTot / system().interaction().sum_inv();
         }
      } else {
         // Explicitly set homogeneous residual components
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (isFlexible()) {
         const int nParam = system().unitCell().nParameter();

         // Combined -1 factor and stress scaling here. This is okay:
         // - residuals only show up as dot products (U, v, norm)
         //   or with their absolute value taken (max), so the
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale,
         //   and then undone right before passing to the unit cell.

         int counter = 0;
         for (int i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               resid[nMonomer*nBasis + counter] = scaleStress_ * -1
                                    * system().mixture().stress(i);
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
      }

   }

   // Update the current system field coordinates
   template <int D>
   void AmIterator<D>::update(DArray<double>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      DArray< DArray<double> > wField;
      wField.allocate(nMonomer);

      // Restructure in format of monomers, basis functions
      for (int i = 0; i < nMonomer; i++) {
         wField[i].allocate(nBasis);
         for (int k = 0; k < nBasis; k++)
         {
            wField[i][k] = newGuess[i*nBasis + k];
         }
      }
      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               wField[i][0] +=
                 system().interaction().chi(i,j) * system().c().basis(j)[0];
            }
         }
      }
      system().setWBasis(wField);

      if (isFlexible()) {
         const int nParam = system().unitCell().nParameter();
         FSArray<double,6> parameters = system().unitCell().parameters();
         int counter = 0;

         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = 1/scaleStress_ * 
                               newGuess[nMonomer*nBasis + counter];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());

         system().setUnitCell(parameters);
      }

   }

   template<int D>
   void AmIterator<D>::outputToLog()
   {
      if (isFlexible()) {
         const int nParam = system().unitCell().nParameter();
         for (int i = 0; i < nParam; i++) {
            Log::file() << "Cell Param  " << i << " = "
                        << Dbl(system().unitCell().parameters()[i], 15)
                        << "  " 
                        << Dbl(system().mixture().stress(i), 15)
                        << "\n";
         }
      }
   }

}
}
#endif
