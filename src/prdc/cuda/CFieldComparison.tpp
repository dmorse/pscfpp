#ifndef PRDC_CUDA_C_FIELD_COMPARISON_TPP
#define PRDC_CUDA_C_FIELD_COMPARISON_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFieldComparison.h"
#include <pscf/cuda/HostDArray.h>
#include <cmath>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   // Default Constructor
   template <int D>
   CFieldComparison<D>::CFieldComparison()
    : maxDiff_(0.0),
      rmsDiff_(0.0)
   {};

   // Comparator for individual fields.
   template <int D>
   double CFieldComparison<D>::compare(CField<D> const& a, 
                                       CField<D> const& b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      int capacity = a.capacity();

      // Allocate arrays on CPU host
      HostDArray<cudaComplex> ha;
      HostDArray<cudaComplex> hb;
      ha.allocate(capacity);
      hb.allocate(capacity);

      // Copy field data from device to host
      ha = a;
      hb = b;

      // Compare fields on host
      double diffSq, diff, d0, d1;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      for (int i = 0; i < capacity; ++i) {
         d0 = ha[i].x - hb[i].x;
         d1 = ha[i].y - hb[i].y;
         diffSq = d0*d0 + d1*d1;
         diff = sqrt(diffSq);
         if (diff > maxDiff_) {
            maxDiff_ = diff;
         }
         rmsDiff_ += diffSq;
      }
      rmsDiff_ = rmsDiff_/double(capacity);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }

   // Comparator for arrays of fields
   template <int D>
   double CFieldComparison<D>::compare(DArray< CField<D> > const & a,
                                       DArray< CField<D> > const & b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      UTIL_CHECK(a[0].capacity() > 0);
      int capacity = a[0].capacity();
      int nFields = a.capacity();

      // Allocate arrays on host
      DArray< HostDArray<cudaComplex> > ha;
      DArray< HostDArray<cudaComplex> > hb;
      ha.allocate(nFields);
      hb.allocate(nFields);
      for (int i = 0; i < nFields; i++) {
         ha[i].allocate(capacity);
         hb[i].allocate(capacity);
      }

      // Copy data from device to host
      for (int i = 0; i < nFields; i++) {
         ha[i] = a[i];
         hb[i] = b[i];
      }

      // Perform comparison using host arrays
      double diffSq, diff, d0, d1;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      int i, j;
      for (i = 0; i < nFields; ++i) {
         for (j = 0; j < capacity; ++j) {
            d0 = ha[i][j].x - hb[i][j].x;
            d1 = ha[i][j].y - hb[i][j].y;
            diffSq = d0*d0 + d1*d1;
            diff = sqrt(diffSq);
            if (diff > maxDiff_) {
               maxDiff_ = diff;
            }
            rmsDiff_ += diffSq;
         }
      }
      rmsDiff_ = rmsDiff_/double(nFields*capacity);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }

}
}
}
#endif
