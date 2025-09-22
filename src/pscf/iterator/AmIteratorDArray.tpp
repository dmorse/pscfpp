#ifndef PSCF_AM_ITERATOR_DARRAY_TPP
#define PSCF_AM_ITERATOR_DARRAY_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NanException.h"
#include <util/global.h> 

namespace Pscf
{

   using namespace Util;

   /*
   * Vector assignment, a = b .
   */
   template <typename Iterator>
   void AmIteratorDArray<Iterator>::setEqual(DArray<double>& a, 
                                             DArray<double> const & b)
   {  a = b; }

   /*
   * Compute and return the inner product of two vectors
   */
   template <typename Iterator> 
   double 
   AmIteratorDArray<Iterator>::dotProduct(DArray<double> const & a,
                                          DArray<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      double product = 0.0;
      for (int i = 0; i < n; i++) {
         // if either value is NaN, throw NanException
         if (std::isnan(a[i]) || std::isnan(b[i])) {
            throw NanException("AmIteratorDArray<Iterator>::dotProduct",
                               __FILE__,__LINE__,0);
         }
         product += a[i] * b[i];
      }
      return product;
   }

   /*
   * Compute and return the maximum magnitude element of a vector.
   */
   template <typename Iterator>
   double AmIteratorDArray<Iterator>::maxAbs(DArray<double> const & a)
   {
      const int n = a.capacity();
      double max = 0.0;
      double value;
      for (int i = 0; i < n; i++) {
         value = a[i];
         if (std::isnan(value)) { // if value is NaN, throw NanException
            throw NanException("AmIteratorDArray<Iterator>::dotProduct",
                                __FILE__,__LINE__,0);
         }
         if (fabs(value) > max)
            max = fabs(value);
      }
      return max;
   }

   /*
   * Compute the vector difference a = b - c .
   */
   template <typename Iterator>
   void AmIteratorDArray<Iterator>::subVV(DArray<double>& a,
                                          DArray<double> const & b,
                                          DArray<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      UTIL_CHECK(n == c.capacity());
      for (int i = 0; i < n; i++) {
         a[i] = b[i] - c[i];
      }
   }

   /*
   * Composite a += b*c for vectors a and b, scalar c .
   */
   template <typename Iterator>
   void AmIteratorDArray<Iterator>::addEqVc(DArray<double>& a,
                                            DArray<double> const & b,
                                            double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      for (int i = 0; i < n; i++) {
         a[i] += c*b[i];
      }
   }

}
#endif
