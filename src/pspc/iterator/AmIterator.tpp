#ifndef PSPC_AM_ITERATOR_TPP
#define PSPC_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "AmIterator.h"
#include <pspc/System.h>
#include <pscf/inter/ChiInteraction.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   template <int D>
   AmIterator<D>::AmIterator(System<D>& system)
   : Iterator<D>(system)
   {}

   template <int D>
   AmIterator<D>::~AmIterator()
   {}

   template <int D>
   double AmIterator<D>::findNorm(FieldCPU const & hist) 
   {
      const int n = hist.capacity();
      double normResSq = 0.0;

      for (int i = 0; i < n; i++) {
         normResSq += hist[i] * hist[i];
      }

      return sqrt(normResSq);
   }

   template <int D>
   double AmIterator<D>::findMaxAbs(FieldCPU const & hist)
   {
      const int n = hist.capacity();
      double maxRes = 0.0;

      for (int i = 0; i < n; i++) {
         if (fabs(hist[i]) > maxRes) 
            maxRes = fabs(hist[i]);
      }

      return maxRes;
   }
   
   template <int D>
   void AmIterator<D>::updateBasis(RingBuffer<FieldCPU> & basis, RingBuffer<FieldCPU> const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      FieldCPU newbasis;
      newbasis.allocate(n);

      for (int i = 0; i < n; i++) {
         newbasis[i] = hists[0][i] - hists[1][i]; // sequential histories basis vectors
      }

      basis.append(newbasis);
   }

   template <int D>
   double AmIterator<D>::computeUDotProd(RingBuffer<FieldCPU> const & resBasis, int m, int n)
   {
      const int length = resBasis[0].capacity();
      
      double dotprod = 0.0;
      for(int i = 0; i < length; i++) {
         dotprod += resBasis[m][i] * resBasis[n][i];
      }

      return dotprod;
   }

   template <int D>
   double AmIterator<D>::computeVDotProd(FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int m)
   {
      const int length = resBasis[0].capacity();
      
      double dotprod = 0.0;
      for(int i = 0; i < length; i++) {
         dotprod += resCurrent[i] * resBasis[m][i];
      }

      return dotprod;
   }

   template <int D>
   void AmIterator<D>::updateU(DMatrix<double> & U, RingBuffer<FieldCPU> const & resBasis, int nHist)
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
   void AmIterator<D>::updateV(DArray<double> & v, FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int nHist)
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }

   template <int D>
   void AmIterator<D>::setEqual(FieldCPU& a, FieldCPU const & b)
   {
      // This seems silly here, but in other implementations it may not be! Goal: no explicit math in AmIteratorTmpl.
      a = b;
   }

   template <int D>
   void AmIterator<D>::addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & basis, DArray<double> coeffs, int nHist)
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
   void AmIterator<D>::addPredictedError(FieldCPU& fieldTrial, FieldCPU const & resTrial, double lambda)
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
   }
   
   template <int D>
   bool AmIterator<D>::hasInitialGuess()
   {
      return sys_->hasWFields();
   }
   
   template <int D>
   int AmIterator<D>::nElements()
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();

      int nEle = nMonomer*nBasis;

      if (sys_->domain().isFlexible()) {
         nEle += sys_->unitCell().nParameter();
      }

      return nEle;
   }

   template <int D>
   void AmIterator<D>::getCurrent(FieldCPU& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();
      const DArray<DArray<double>> * currSys = &sys_->wFields(); 
      
      
      for (int i = 0; i < nMonomer; i++) {
         for (int k = 0; k < nBasis; k++)
         {
            curr[i*nBasis+k] = (*currSys)[i][k];
         }
      }

      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         const FSArray<double,6> currParam = sys_->unitCell().parameters();

         for (int i = 0; i < nParam; i++) {
            curr[nMonomer*nBasis + i] = scaleStress*currParam[i];
         }
      }

      return;
   }

   template <int D>
   void AmIterator<D>::evaluate()
   {
      // Solve MDEs for current omega field
      sys_->compute();
      // Compute stress if done
      if (sys_->domain().isFlexible()) {
         sys_->mixture().computeStress();
      }
   }

   template <int D>
   void AmIterator<D>::getResidual(FieldCPU& resid)
   {
      const int n = nElements();
      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();

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
                  sys_->interaction().chi(i,j)*sys_->cField(j)[k] -
                  sys_->interaction().idemp(i,j)*sys_->wField(j)[k];
            }
         }
      }

      // If not canonical, account for incompressibility 
      if (!sys_->mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= 1.0/sys_->interaction().sum_inv();
         }
      } else {
         // otherwise explicitly set the residual value for the homogeneous components
         // the homogeneous field component is set explicitly as well in update 
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         
         // Combined -1 factor and stress scaling here. This is okay: 
         // - residuals only show up as dot products (U, v, norm) 
         //   or with their absolute value taken (max), so the 
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale, 
         //   and then undone right before passing to the unit cell.

         for (int i = 0; i < nParam ; i++) {
            resid[nMonomer*nBasis + i] = scaleStress * -1 
                                       * sys_->mixture().stress(i);
         }
      }

   }

   template <int D>
   void AmIterator<D>::update(FieldCPU& newGuess)
   {
      // Convert back to field format
      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();
      
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
      // Manually and explicitly set homogeneous components of field if canonical
      if (sys_->mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               wField[i][0] += 
                  sys_->interaction().chi(i,j) * sys_->cField(j)[0];
            }
         }
      }
      sys_->setWBasis(wField);

      if (sys_->domain().isFlexible()) {
         FSArray<double, 6> parameters;
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();

         for (int i = 0; i < nParam; i++) {
            parameters.append(1/scaleStress * newGuess[nMonomer*nBasis + i]);
         }

         sys_->setUnitCell(parameters);
      }
      
   }

   template<int D>
   void AmIterator<D>::outputToLog()
   {
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         for (int i = 0; i < nParam; i++) {
            Log::file() << "Parameter " << i << " = "
                        << Dbl(sys_->unitCell().parameters()[i])
                        << "\n";
         }
      }
   }



}
}
#endif
