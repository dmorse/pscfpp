#ifndef PSPG_AM_ITERATOR_BASIS_TPP
#define PSPG_AM_ITERATOR_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include <pspg/System.h>
#include <pspg/field/RDField.h>

#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>

#include <pscf/inter/Interaction.h>
#include <pscf/iterator/NanException.h>

#include <util/global.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;
   using namespace Prdc;

   // Constructor
   template <int D>
   AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
   : Iterator<D>(system)
   {  setClassName("AmIteratorBasis"); }

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

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);
   }

   // -- Protected virtual function -- //

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorBasis<D>::setup(bool isContinuation)
   {
      // Setup by AM algorithm
      AmIteratorTmpl<Iterator<D>, DArray<double> >::setup(isContinuation);

      // Update chi matrix and related properties in member interaction_
      interaction_.update(system().interaction());
   }


   // -- Private virtual functions used to implement AM algorithm --  //

   template <int D>
   void AmIteratorBasis<D>::setEqual(DArray<double>& a, 
                                     DArray<double> const & b)
   {  a = b; }

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

   #if 0
   template <int D>
   double AmIteratorBasis<D>::norm(DArray<double>  const & a) 
   {
      const int n = a.capacity();
      double data;
      double normSq = 0.0;
      for (int i=0; i < n; ++i) {
         data = a[i];
         normSq += data*data;
      }
      return sqrt(normSq);
   }

   // Compute one element of U matrix of by computing a dot product
   template <int D>
   double
   AmIteratorBasis<D>::computeUDotProd(RingBuffer<DArray<double> > const & resBasis,
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
   AmIteratorBasis<D>::computeVDotProd(DArray<double> const & resCurrent,
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
   void AmIteratorBasis<D>::updateU(DMatrix<double> & U,
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
   void AmIteratorBasis<D>::updateV(DArray<double> & v,
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
   void 
   AmIteratorBasis<D>::updateBasis(RingBuffer<DArray<double> > & basis,
                                   RingBuffer<DArray<double> > const& hists)
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

   // Compute and return the number of elements in a field vector
   template <int D>
   int AmIteratorBasis<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      int nEle = nMonomer*nBasis;

      if (isFlexible_) { 
         nEle += system().unitCell().nParameter();
      }

      return nEle;
   }

   // Get the current field from the system
   template <int D>
   void AmIteratorBasis<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      const DArray< DArray<double> > * currSys = &system().w().basis();

      for (int i = 0; i < nMonomer; i++) {
         for (int k = 0; k < nBasis; k++) {
            curr[i*nBasis+k] = (*currSys)[i][k];
         }
      }

      if (isFlexible_) {
         const int begin = nMonomer*nBasis;
         const int nParam = system().unitCell().nParameter();
         FSArray<double,6> const & parameters
                                  = system().unitCell().parameters();
         for (int i = 0; i < nParam; i++) {
            curr[begin + i] = scaleStress_ * parameters[i];
         }
      }

   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void AmIteratorBasis<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute();

      // If flexible, compute stress 
      if (isFlexible_) {
         system().mixture().computeStress(system().domain().waveList());
      }
   }

   // Compute the residual for the current system state
   template <int D>
   void AmIteratorBasis<D>::getResidual(DArray<double>& resid)
   {
      UTIL_CHECK(system().basis().isInitialized());
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      const int n = nElements();

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
                  interaction_.chi(i,j)*system().c().basis(j)[k] -
                  interaction_.p(i,j)*system().w().basis(j)[k];
            }
         }
      }

      // If not canonical, account for incompressibility
      if (!system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= 1.0/interaction_.sumChiInverse();
         }
      } else {
         // Explicitly set homogeneous residual components
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         for (int i = 0; i < nParam ; i++) {
            resid[nMonomer*nBasis + i] = -1.0 * scaleStress_ 
                                       * system().mixture().stress(i);
         }

         //  Note: 
         //  Combined -1 factor and stress scaling here.  This is okay:
         //  - Residuals only show up as dot products (U, v, norm)
         //    or with their absolute value taken (max), so the
         //    sign on a given residual vector element is not relevant
         //    as long as it is consistent across all vectors
         //  - The scaling is applied here and to the unit cell param
         //    storage, so that updating is done on the same scale,
         //    and then undone right before passing to the unit cell.
      }

   }

   // Update the current system field coordinates
   template <int D>
   void AmIteratorBasis<D>::update(DArray<double>& newGuess)
   {
      UTIL_CHECK(system().basis().isInitialized());
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

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
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               wField[i][0] +=
                 interaction_.chi(i,j) * system().c().basis(j)[0];
            }
         }
      }
      system().setWBasis(wField);

      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         const int begin = nMonomer*nBasis;
         FSArray<double,6> parameters;
         double value;
         for (int i = 0; i < nParam; i++) {
            value = newGuess[begin + i]/scaleStress_;
            parameters.append(value);
         }
         system().setUnitCell(parameters);
      }

   }

   template<int D>
   void AmIteratorBasis<D>::outputToLog()
   {
      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         for (int i = 0; i < nParam; i++) {
            Log::file() << "Parameter " << i << " = "
                        << Dbl(system().unitCell().parameters()[i])
                        << "\n";
         }
      }
   }

}
}
#endif
