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

namespace Pscf {
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
   {  setClassName("AmIterator"); }

   // Read parameters from file
   template <int D>
   void AmIterator<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Iterator<D>, DArray<double> >::readParameters(in);

      // Default parameter values
      bool isFlexible = 0; 
      scaleStress_ = 10.0;

      // Setup to read optional flexibleParams boolean array
      DArray<int> flexParamBools;
      int np = system().unitCell().nParameter();
      flexParamBools.allocate(np);
      flexParamBools[0] = -1; // Flag to determine if array has been set

      // Read in additional optional parameters
      readOptional(in, "isFlexible", isFlexible);
      readOptionalDArray(in, "flexibleParams", flexParamBools, np);
      readOptional(in, "scaleStress", scaleStress_);

      // Set flexibleParams_ array based on user input
      if (flexParamBools[0] != -1) { 
         // user declared flexibleParams explicitly
         for (int i = 0; i < np; i++) {
            if (flexParamBools[i]) {
               flexibleParams_.append(i);
            }
         }
      } else if (isFlexible) {
         // all parameters are flexible
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(i);
         }
      }
   }

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

   // Compute and return maximum element of residual vector.
   template <int D>
   double AmIterator<D>::findMaxAbs(DArray<double> const & hist)
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

   template <int D>
   void AmIterator<D>::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

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

   // Does the system have an initial field guess?
   template <int D>
   bool AmIterator<D>::hasInitialGuess()
   {
      return system().w().hasData();
   }

   // Compute and return the number of elements in a field vector
   template <int D>
   int AmIterator<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      int nEle = nMonomer*nBasis;

      nEle += flexibleParams().size();

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

      const FSArray<int,6> indices = flexibleParams();
      const int nParam = indices.size();
      const FSArray<double,6> currParam = system().unitCell().parameters();

      for (int i = 0; i < nParam; i++) {
         curr[nMonomer*nBasis + i] = scaleStress_*currParam[indices[i]];
      }

      return;
   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void AmIterator<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute();
      // Compute stress if done
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
         const FSArray<int,6> indices = flexibleParams();
         const int nParam = indices.size();

         // Combined -1 factor and stress scaling here. This is okay:
         // - residuals only show up as dot products (U, v, norm)
         //   or with their absolute value taken (max), so the
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale,
         //   and then undone right before passing to the unit cell.

         for (int i = 0; i < nParam ; i++) {
            resid[nMonomer*nBasis + i] = scaleStress_ * -1
                                       * system().mixture().stress(indices[i]);
         }
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
         const FSArray<int,6> indices = flexibleParams();
         const int nParam = indices.size();
         FSArray<double,6> parameters = system().unitCell().parameters();
         int ind, i;

         for (i = 0; i < nParam; i++) {
            ind = indices[i];
            parameters[ind] = 1/scaleStress_ * newGuess[nMonomer*nBasis + i];
         }

         system().setUnitCell(parameters);
      }

   }

   template<int D>
   void AmIterator<D>::outputToLog()
   {
      if (isFlexible()) {
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
