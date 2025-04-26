#ifndef PRDC_CUDA_R_FIELD_DFT_COMPARISON_TPP
#define PRDC_CUDA_R_FIELD_DFT_COMPARISON_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDftComparison.h"
#include <pscf/cuda/HostDArray.h>
#include <cmath>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   // Default Constructor
   template <int D>
   RFieldDftComparison<D>::RFieldDftComparison()
    : maxDiff_(0.0),
      rmsDiff_(0.0)
   {};

   // Comparator for individual fields.
   template <int D>
   double RFieldDftComparison<D>::compare(RFieldDft<D> const& a, 
                                          RFieldDft<D> const& b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      int capacity = a.capacity();

      HostDArray<cudaComplex> ha;
      HostDArray<cudaComplex> hb;
      ha.allocate(capacity);
      hb.allocate(capacity);
      ha = a;
      hb = b;

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
   double RFieldDftComparison<D>::compare(DArray< RFieldDft<D> > const & a,
                                          DArray< RFieldDft<D> > const & b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      UTIL_CHECK(a[0].capacity() > 0);
      int capacity = a[0].capacity();
      int nFields = a.capacity();

      DArray< HostDArray<cudaComplex> > ha;
      DArray< HostDArray<cudaComplex> > hb;
      ha.allocate(nFields);
      hb.allocate(nFields);
      for (int i = 0; i < nFields; i++) {
         ha[i].allocate(capacity);
         hb[i].allocate(capacity);
         ha[i] = a[i];
         hb[i] = b[i];
      }

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
