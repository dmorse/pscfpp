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
#include <pscf/inter/Interaction.h>
#include <pspg/field/RDField.h>
#include <util/global.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;

   template <int D>
   AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
   : Iterator<D>(system)
   { setClassName("AmIteratorBasis"); }

   template <int D>
   AmIteratorBasis<D>::~AmIteratorBasis()
   {}

   template <int D>
   void AmIteratorBasis<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters 
      AmIteratorTmpl<Iterator<D>,DArray<double> >::readParameters(in);

      // Default parameter values
      isFlexible_ = 0;
      scaleStress_ = 10.0;

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);
   }


   template <int D>
   double AmIteratorBasis<D>::findNorm(DArray<double>  const & hist) 
   {
      const int n = hist.capacity();
      double data;
      double normResSq = 0.0;
      for (int i=0; i < n; ++i) {
         data = hist[i];
         normResSq += data*data;
      } 
      return sqrt(normResSq);
   }

   // Compute and return maximum element of residual vector.
   template <int D>
   double AmIteratorBasis<D>::findMaxAbs(DArray<double> const & hist)
   {
      const int n = hist.capacity();
      double maxRes = 0.0;

      for (int i = 0; i < n; i++) {
         if (fabs(hist[i]) > maxRes)
            maxRes = fabs(hist[i]);
      }

      return maxRes;
   }

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

   template <int D>
   void AmIteratorBasis<D>::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

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
         system().mixture().computeStress(system().wavelist());
      }
   }

   // Compute the residual for the current system state
   template <int D>
   void AmIteratorBasis<D>::getResidual(DArray<double>& resid)
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

      // If not canonical, account for incompressibility
      if (!system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= 1.0/system().interaction().sum_inv();
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
