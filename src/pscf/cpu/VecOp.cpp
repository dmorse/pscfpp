/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "VecOp.h"
#include <util/containers/Array.h>
#include <cmath>

namespace Pscf {
namespace VecOp {

   // Assignment

   /*
   * Vector-vector assignment, a[i] = b[i] (real, slice).
   */
   void eqV(Array<double>& a, Array<double> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(beginIdA >= 0);
      UTIL_CHECK(beginIdB >= 0);
      UTIL_CHECK(n > 0);
      UTIL_CHECK(a.capacity() >= beginIdA + n);
      UTIL_CHECK(b.capacity() >= beginIdB + n);
      for (int i = 0; i < n; ++i) {
         a[i + beginIdA] = b[i + beginIdB];
      }
   }

   /*
   * Vector-vector assignment, a[i] = b[i] (real).
   */
   void eqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i];
      }
   }

   /*
   * Vector-scalar assignment, a[i] = b (real).
   */
   void eqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] = b;
      }
   }

   // Addition

   /*
   * Vector-vector addition, a[i] = b[i] + c[i] (real).
   */
   void addVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] + c[i];
      }
   }

   /*
   * Vector-scalar addition, a[i] = b[i] + c (real).
   */
   void addVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);

      for (int i = 0; i < n; ++i) {
         a[i] = b[i] + c;
      }
   }

   // Subtraction

   /*
   * Vector-vector subtraction, a[i] = b[i] - c[i] (real).
   */
   void subVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] - c[i];
      }
   }

   /*
   * Vector-scalar subtraction, a[i] = b[i] - c (real).
   */
   void subVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] - c;
      }
   }

   // Multiplication

   /*
   * Vector-vector multiplication, a[i] = b[i] * c[i] (real).
   */
   void mulVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] * c[i];
      }
   }

   /*
   * Vector-scalar multiplication, a[i] = b[i] * c (real).
   */
   void mulVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] * c;
      }
   }

   // Division

   /*
   * Vector-vector division, a[i] = b[i] / c[i] (real).
   */
   void divVV(Array<double>& a,
              Array<double> const & b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] / c[i];
      }
   }

   /*
   * Vector-scalar division, a[i] = b[i] / c (real).
   */
   void divVS(Array<double>& a, Array<double> const & b, double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i] / c;
      }
   }

   /*
   * Scalar-vector division, a[i] = b / c[i] (real).
   */
   void divSV(Array<double>& a, double b, Array<double> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(c.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b/c[i];
      }
   }

   // In-place addition

   /*
   * Vector-vector in-place addition, a[i] += b[i] (real).
   */
   void addEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] += b[i];
      }
   }

   /*
   * Vector-scalar in-place addition, a[i] += b (real).
   */
   void addEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] += b;
      }
   }

   // In-place subtraction

   /*
   * Vector-vector in-place subtraction, a[i] -= b[i] (real).
   */
   void subEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] -= b[i];
      }
   }

   /*
   * Vector-scalar in-place subtraction, a[i] -= b (real).
   */
   void subEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] -= b;
      }
   }

   // In-place multiplication

   /*
   * Vector-vector in-place multiplication, a[i] *= b[i] (real).
   */
   void mulEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] *= b[i];
      }
   }

   /*
   * Vector-scalar in-place multiplication, a[i] *= b (real).
   */
   void mulEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] *= b;
      }
   }

   // In-place division

   /*
   * Vector-vector in-place division, a[i] /= b[i] (real).
   */
   void divEqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] /= b[i];
      }
   }

   /*
   * Vector-scalar in-place division, a[i] /= b (real).
   */
   void divEqS(Array<double>& a, double b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      for (int i = 0; i < n; ++i) {
         a[i] /= b;
      }
   }

   // Exponentiation

   /*
   * Vector exponentiation, a[i] = exp(b[i]) (real).
   */
   void expV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = exp(b[i]);
      }
   }

   /*
   * Exponentiate a scaled vector, a[i] = exp(b[i]*c) (real).
   */
   void expVc(Array<double>& a, Array<double> const & b, const double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = exp(b[i]*c);
      }
   }

   // Square

   /*
   * Vector elementwise square, a[i] = b[i]*b[i] (real).
   */
   void sqV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i]*b[i];
      }
   }

   // Absolute magnitude

   /*
   * Vector elementwise square, a[i] = b[i]*b[i] (real).
   */
   void absV(Array<double>& a, Array<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = std::abs(b[i]);
      }
   }

   // Linear combinations (separate result)

   /*
   * Add two scaled vectors, a[i] = b1[i]*c1 + b2[2]*c2 (real).
   */
   void addVcVc(Array<double>& a,
                Array<double> const & b1, const double c1,
                Array<double> const & b2, const double c2)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b1[i]*c1 + b2[i]*c2;
      }
   }

   /*
   * Add a scaled vector and a scalar, a[i] = b[i]*c + s (real).
   */
   void addVcS(Array<double>& a,
               Array<double> const & b, const double c,
               const double s)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b[i]*c + s;
      }
   }

   /*
   * In-place add one scaled vector, a[i] += b[i]*c (real).
   */
   void addEqVc(Array<double>& a, 
                Array<double> const & b, const double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] += b[i]*c;
      }
   }

   /*
   * Add 2 scaled vectors & scalar, a[i] = b1[i]*c1 + b2[i]*c2 + s (real).
   */
   void addVcVcS(Array<double>& a,
                 Array<double> const & b1, const double c1,
                 Array<double> const & b2, const double c2,
                 const double s)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b1[i]*c1 + b2[i]*c2 + s;
      }
   }

   /*
   * Add 3 scaled vectors, a[i] = b1[i]*c1 + b2[i]*c2 + b3[i]*c3 (real).
   */
   void addVcVcVc(Array<double>& a,
                  Array<double> const & b1, const double c1,
                  Array<double> const & b2, const double c2,
                  Array<double> const & b3, const double c3)
   {
      const int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);
      UTIL_CHECK(b3.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a[i] = b1[i]*c1 + b2[i]*c2 + b3[i]*c3;
      }
   }

   // Pair operations (two output arrays and a shared input)

   /*
   * Vector assignment in pairs, ax[i] = b[i], x = 1, 2 (real).
   */
   void eqVPair(Array<double>& a1, 
                Array<double>& a2,
                Array<double> const & b)
   {
      const int n = a1.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(a2.capacity() == n);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a1[i] = b[i];
         a2[i] = b[i];
      }
   }

   /*
   * Vector multiplication in pairs, ax[i] = bx[i] * c[i], x=1,2 (real).
   */
   void mulVVPair(Array<double>& a1, 
                  Array<double>& a2,
                  Array<double> const & b1,
                  Array<double> const & b2,
                  Array<double> const & c)
   {
      const int n = a1.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(a2.capacity() == n);
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a1[i] = b1[i]*c[i];
         a2[i] = b2[i]*c[i];
      }
   }

   /*
   * In-place vector multiplication in pairs, ax[i] *= b[i], x=1,2 (real).
   */
   void mulEqVPair(Array<double>& a1, 
                   Array<double>& a2,
                   Array<double> const & b)
   {
      const int n = a1.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(a2.capacity() == n);
      UTIL_CHECK(b.capacity() >= n);
      for (int i = 0; i < n; ++i) {
         a1[i] *= b[i];
         a2[i] *= b[i];
      }
   }

} // namespace VecOp
} // namespace Pscf
