#ifndef PRDC_CUDA_K_FIELD_COMPARISON_TPP
#define PRDC_CUDA_K_FIELD_COMPARISON_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "KFieldComparison.h"
#include <cmath>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   // Default Constructor
   template <int D>
   KFieldComparison<D>::KFieldComparison()
    : maxDiff_(0.0),
      rmsDiff_(0.0)
   {};

   // Comparator for individual fields.
   template <int D>
   double KFieldComparison<D>::compare(RFieldDft<D> const& a, RFieldDft<D> const& b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());

      // Create temporary host arrays
      int nPoints = a.capacity();
      cudaComplex* temp_a = new cudaComplex[nPoints];
      cudaComplex* temp_b = new cudaComplex[nPoints];
      cudaMemcpy(temp_a, a.cField(), nPoints*sizeof(cudaComplex), cudaMemcpyDeviceToHost);
      cudaMemcpy(temp_b, b.cField(), nPoints*sizeof(cudaComplex), cudaMemcpyDeviceToHost);

      double diffSq, diff, d0, d1;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      for (int i = 0; i < nPoints; ++i) {
         d0 = temp_a[i].x - temp_b[i].x;
         d1 = temp_a[i].y - temp_b[i].y;
         diffSq = d0*d0 + d1*d1;
         diff = sqrt(diffSq);
         if (diff > maxDiff_) {
            maxDiff_ = diff;
         }
         rmsDiff_ += diffSq;
      }
      rmsDiff_ = rmsDiff_/double(nPoints);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }

   // Comparator for arrays of fields
   template <int D>
   double KFieldComparison<D>::compare(DArray< RFieldDft<D> > const & a,
                                       DArray< RFieldDft<D> > const & b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      UTIL_CHECK(a[0].capacity() > 0);

      // Create temporary host arrays
      DArray< cudaComplex* > temp_a;
      DArray< cudaComplex* > temp_b;

      int nFields = a.capacity();
      int nPoints = a[0].capacity();
      temp_a.allocate(nFields);
      temp_b.allocate(nFields);
      for (int i = 0; i < nFields; i++) {
         temp_a[i] = new cudaComplex[nPoints];
         temp_b[i] = new cudaComplex[nPoints];
         cudaMemcpy(temp_a[i], a[i].cField(), nPoints*sizeof(cudaComplex), cudaMemcpyDeviceToHost);
         cudaMemcpy(temp_b[i], b[i].cField(), nPoints*sizeof(cudaComplex), cudaMemcpyDeviceToHost);
      }      

      double diffSq, diff, d0, d1;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      int i, j;
      for (i = 0; i < nFields; ++i) {
         for (j = 0; j < nPoints; ++j) {
            d0 = temp_a[i][j].x - temp_b[i][j].x;
            d1 = temp_a[i][j].y - temp_b[i][j].y;
            diffSq = d0*d0 + d1*d1;
            diff = sqrt(diffSq);
            if (diff > maxDiff_) {
               maxDiff_ = diff;
            }
            rmsDiff_ += diffSq;
         }
      }
      rmsDiff_ = rmsDiff_/double(nFields*nPoints);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }

}
}
}
#endif
