#ifndef PSCF_FIELD_COMPARISON_TPP
#define PSCF_FIELD_COMPARISON_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldComparison.h"
#include <cmath>

namespace Pscf {

   using namespace Util;

   // Default Constructor
   template <class FT>
   FieldComparison<FT>::FieldComparison(int begin)
    : maxDiff_(0.0),
      rmsDiff_(0.0),
      begin_(begin)
   {};

   // Comparator for individual fields.
   template <class FT>
   double FieldComparison<FT>::compare(FT const& a, FT const& b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      int n = a.capacity();
      double diff;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      for (int i = begin_; i < n; ++i) {
         diff = std::abs(a[i] - b[i]);
         if (std::isnan(diff)) {
            // If either field has a NaN component, set error to very
            // high value and exit the function
            maxDiff_ = 1e8;
            rmsDiff_ = 1e8;
            return maxDiff_;
         } else if (diff > maxDiff_) {
            maxDiff_ = diff;
         }
         rmsDiff_ += diff*diff;
      }
      rmsDiff_ = rmsDiff_/double(n);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }

   // Comparator for arrays of fields
   template <class FT>
   double FieldComparison<FT>::compare(DArray<FT> const & a,
                                       DArray<FT> const & b)
   {
      UTIL_CHECK(a.capacity() > 0);
      UTIL_CHECK(a.capacity() == b.capacity());
      UTIL_CHECK(a[0].capacity() > 0);
      int m = a.capacity();
      double diff;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      int i, j, n;
      for (i = 0; i < m; ++i) {
         n = a[i].capacity();
         UTIL_CHECK(n > 0);
         UTIL_CHECK(n == b[i].capacity());
         for (j = begin_; j < n; ++j) {
            diff = std::abs(a[i][j] - b[i][j]);
            if (std::isnan(diff)) {
               // If either field has a NaN component, set error to very
               // high value and exit the function
               maxDiff_ = 1e8;
               rmsDiff_ = 1e8;
               return maxDiff_;
            } else if (diff > maxDiff_) {
               maxDiff_ = diff;
            }
            rmsDiff_ += diff*diff;
         }
      }
      rmsDiff_ = rmsDiff_/double(m*n);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }

}
#endif
