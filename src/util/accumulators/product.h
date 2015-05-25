#ifndef UTIL_PRODUCT_H
#define UTIL_PRODUCT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <complex>
using std::complex;

namespace Util
{

   /**
   * Product for float Data.
   */
   inline float product(float a, float b)
   {  return a*b; }

   /**
   * Product for double Data.
   */
   inline double product(double a, double b)
   {  return a*b; }

   /**
   * Dot product for Vector Data.
   */
   inline double product(const Vector& a, const Vector& b)
   {  return a.dot(b); }

   /**
   * Double contraction for Tensor Data.
   */
   inline double product(const Tensor& a, const Tensor& b)
   {
      double sum = 0.0;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            sum += a(i, j)*b(i, j);
         }
      }
      return sum;
   }

   /**
   * Inner product for complex<float> Data.
   */
   inline complex<float> product(complex<float> a, complex<float> b)
   {  return conj(a)*b; }

   /**
   * Inner product for complex<double> Data.
   */
   inline complex<double> product(complex<double> a, complex<double> b)
   {  return conj(a)*b; }

}
#endif
