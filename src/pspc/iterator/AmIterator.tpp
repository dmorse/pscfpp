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
   {  setClassName("AmIterator"); }

   template <int D>
   AmIterator<D>::~AmIterator()
   {  setClassName("AmIterator"); }

   template <int D>
   void AmIterator<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters 
      AmIteratorTmpl<Iterator<D>,FieldCPU>::readParameters(in);

      // Default parameter values
      isFlexible_ = 0;
      scaleStress_ = 10.0;

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);

      // If lattice parameters are flexible, update flexibleParams_
      if (isFlexible_) {
         // all parameters are flexible
         int np = system().unitCell().nParameter();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(i);
         }
      }
   }

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
      return system().hasWFields();
   }
   
   template <int D>
   int AmIterator<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      int nEle = nMonomer*nBasis;

      nEle += flexibleParams().size();

      return nEle;
   }

   template <int D>
   void AmIterator<D>::getCurrent(FieldCPU& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      const DArray<FieldCPU> * currSys = &system().wFields(); 
      
      
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

   template <int D>
   void AmIterator<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute();
      // Compute stress if done
      if (isFlexible_) {
         system().mixture().computeStress();
      }
   }

   template <int D>
   void AmIterator<D>::getResidual(FieldCPU& resid)
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
                  system().interaction().chi(i,j)*system().cField(j)[k] -
                  system().interaction().idemp(i,j)*system().wField(j)[k];
            }
         }
      }

      // If iterator has a mask, account for the mask in the values of the residuals
      if (hasMask()) {
         for (int i = 0; i < nMonomer; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] += maskField()[k] / system().interaction().sum_inv(); 
            }
         }
      }

      // If not canonical, account for incompressibility 
      if (!system().mixture().isCanonical()) {
         double phi; // volume fraction of mask, if a mask exists
         if (hasMask()) {
            phi = maskField()[0];
         } else {
            phi = 0.0;
         }

         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= (1.0-phi)/system().interaction().sum_inv();
         }
      } else {
         // otherwise explicitly set the residual value for the homogeneous components
         // the homogeneous field component is set explicitly as well in update 
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (isFlexible_) {
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

   template <int D>
   void AmIterator<D>::update(FieldCPU& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      
      DArray<FieldCPU> wField;
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
      if (system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               wField[i][0] += 
                  system().interaction().chi(i,j) * system().cField(j)[0];
            }
         }
      }
      system().setWBasis(wField);

      if (isFlexible_) {
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
