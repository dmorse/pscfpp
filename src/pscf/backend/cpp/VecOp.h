#ifndef PSCF_CPU_VEC_OP_H
#define PSCF_CPU_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Element-wise vector operations on Cpu arrays.
   *
   * Operations that are performed by these functions include addition,
   * subtraction, multiplication, division, exponentiation, and assignment.
   * The function names will, correspondingly, begin with "add", "sub",
   * "mul", "div", "exp", or "eq" to indicate the relevant operation.
   * Functions are also included to perform compound or in-place arithmetic
   * operations, i.e.  those that are performed using +=, -=, *=, and /= in
   * C++. These functions have names that begin with "addEq", "subEq",
   * "mulEq", and "divEq", respectively.
   *
   * The output (the LHS of the vector operation) is always the first
   * parameter passed to the function. The input argument(s) (on the RHS
   * of the vector operation) may be vectors or scalars. If an argument
   * is a vector (scalar), the function name will contain a V (S). The
   * order in which input parameters are listed in a function interface
   * is always the same as the order in which the corresponding symbols
   * V and/or S appear in the function name. For example, the function
   * addVV(A,B,C) implements vector-vector addition A[i] = B[i] + C[i],
   * while addVS(A,B,c) implements vector-scalar addition A[i] = B[i] + c
   * in which c is a scalar that is added to every element of B. For
   * commutative operations involving both vectors and scalars, the
   * vectors are listed first. So, for example, addVS exists, but addSV
   * does not.
   *
   * Some functions use the product of a vector and a scalar coefficient
   * as an inut to a subsequent operation, such as addition of two or more
   * such "scaled" vectors or exponentiation of a scaled vector. In the
   * names of such functions, the symbol "Vc" is used to indicate an input
   * vector that is multiplied by a corresponding scalar coefficient. When
   * such a vector-scalar product is indicated, the input vector and
   * corresponding scalar coefficient always appear next to each other as
   * parameters in the function interface, with the vector listed before
   * the scalar coefficient.  For example, the function
   * addVcVc(a, b1, c1, b2, c2), in which a, b1 and b2 are arrays while
   * c1 and c2 are scalar coefficients, performs the linear combination
   * a[i] = b1[i]*c1 + b2[i}*c2 for every value of the element index i.
   *
   * \defgroup Pscf_Cpu_VecOp_Module VecOp (CPU)
   * \ingroup Pscf_Cpu_Module
   */
   namespace VecOp {

      /*
      * This file and VecOp.cpp only contain declarations and definitions
      * for operations that involve only real arrays and scalars.
      * Corresponding declarations and definitions for vector operations
      * that involve complex-valued arrays and/or scalars are given in
      * the files VecOpCx.h and VecOpCx.cpp, respectively.
      */

      // Assignment

      /**
      * Vector assignment, a[i] = b[i] (real, slice).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in slice of array a (HHS)
      * \param beginIdB  index of first element in slice of array b (RHS)
      * \param n  number of elements in the slice
      */
      void eqV(Array<double>& a, Array<double> const & b,
               const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector assignment, a[i] = b[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void eqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector assignment, a[i] = b (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      void eqS(Array<double>& a, double b);

      // Addition

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (real)
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      void addVV(Array<double>& a,
                 Array<double> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      void addVS(Array<double>& a, Array<double> const & b, double c);

      // Subtraction

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (real)
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      void subVV(Array<double>& a,
                 Array<double> const & b, Array<double> const & c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      void subVS(Array<double>& a, Array<double> const & b, double c);

      // Multiplication

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      void mulVV(Array<double>& a,
                 Array<double> const & b, Array<double> const & c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      void mulVS(Array<double>& a, Array<double> const & b, double c);

      // Division

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      void divVV(Array<double>& a,
                 Array<double> const & b, Array<double> const & c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      void divVS(Array<double>& a, Array<double> const & b, double c);

      /**
      * Vector division, a[i] = b / c[i].
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param c  real array (RHS)
      */
      void divSV(Array<double>& a, double b, Array<double> const & c);

      // In-place addition

      /**
      * Vector-vector in-place addition, a[i] += b[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void addEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place addition, a[i] += b (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      void addEqS(Array<double>& a, double b);

      // In-place subtraction

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void subEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar subtraction in-place, a[i] -= b (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      void subEqS(Array<double>& a, double b);

      // In-place multiplication

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void mulEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      void mulEqS(Array<double>& a, double b);

      // In-place division

      /**
      * Vector-vector in-place division, a[i] /= b[i].
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void divEqV(Array<double>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place division, a[i] /= b.
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      void divEqS(Array<double>& a, double b);

      // Exponentiation

      /**
      * Vector exponentiation, a[i] = exp(b[i]) (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void expV(Array<double>& a, Array<double> const & b);

      /**
      * Exponentiation a scaled vector, a[i] = exp(b[i]*c) (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar coefficient of b (RHS)
      */
      void expVc(Array<double>& a,
                 Array<double> const & b, const double c);

      // Square

      /**
      * Vector element-wise square, a[i] = b[i]*b[i] (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void sqV(Array<double>& a, Array<double> const & b);

      // Absolute magnitude

      /**
      * Element-wise absolute magnitude, a[i] = abs(b[i]) (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      void absV(Array<double>& a, Array<double> const & b);

      // Addition of scaled vectors (linear combinations)

      /**
      * Add two scaled vectors, a[i] = b1[i]*c1 + b2[2]*c2 (real).
      *
      * \param a  real output array (LHS)
      * \param b1  1st real input array (RHS)
      * \param c1  real coefficient of b1 (RHS)
      * \param b2  2nd real input array (RHS)
      * \param c2  real coefficient of b2 (RHS)
      */
      void addVcVc(Array<double>& a,
                   Array<double> const & b1, const double c1,
                   Array<double> const & b2, const double c2);

      /**
      * Add a scaled vector and a scalar, a[i] = b[i]*c + s (real).
      *
      * \param a  real output array (LHS)
      * \param b  real input array (RHS)
      * \param c  real coefficient of b (RHS)
      * \param s  real scalar summand (RHS)
      */
      void addVcS(Array<double>& a,
                  Array<double> const & b, const double c,
                  const double s);

      /**
      * Add scaled vector in-place, a[i] += b[i]*c (real).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real coefficient of b (RHS)
      */
      void addEqVc(Array<double>& a,
                   Array<double> const & b, double const c);

      /**
      * Add scaled vectors + scalar, a[i] = b1[i]*c1 + b2[2]*c2 + s (real).
      *
      * \param a  real array (LHS)
      * \param b1  1st real input array (RHS)
      * \param c1  real coefficient of b1 (RHS)
      * \param b2  2nd real input array (RHS)
      * \param c2  real coefficient of b2 (RHS)
      * \param s  real scalar summand (RHS)
      */
      void addVcVcS(Array<double>& a,
                   Array<double> const & b1, const double c1,
                   Array<double> const & b2, const double c2,
                   const double s);

      /**
      * Add scaled vectors, a[i] = b1[i]*c1 + b2[i]*c2 + b3[i]*c3 (real).
      *
      * \param a  real array (LHS)
      * \param b1  1st real input array (RHS)
      * \param c1  real coefficient of b1 (RHS)
      * \param b2  2nd real input array (RHS)
      * \param c2  real coefficient of b2 (RHS)
      * \param b3  3rd real input array (RHS)
      * \param c3  real coefficient of b3 (RHS)
      */
      void addVcVcVc(Array<double>& a,
                     Array<double> const & b1, const double c1,
                     Array<double> const & b2, const double c2,
                     Array<double> const & b3, const double c3);

      // Pair operations (two output arrays and a shared input)

      /**
      * Vector assignment in pairs, ax[i] = b[i], x = 1, 2.
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a1  real array 1 (LHS)
      * \param a2  real array 2 (LHS)
      * \param b  shared real array to be assigned to both a1 and a2
      */
      void eqVPair(Array<double>& a1,
                   Array<double>& a2,
                   Array<double> const & b);

      /**
      * Vector multiplication in pairs, ax[i] = bx[i] * s[i], x=1,2.
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a1  real array 1 (LHS)
      * \param a2  real array 2 (LHS)
      * \param b1  real array 1 (RHS)
      * \param b2  real array 2 (RHS)
      * \param c  shared real array to multiply both b1 and b2
      */
      void mulVVPair(Array<double>& a1, Array<double>& a2,
                     Array<double> const & b1,
                     Array<double> const & b2,
                     Array<double> const & c);

      /**
      * In-place vector multiplication in pairs, ax[i] *= b[i], x=1,2.
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a1  real array 1 (LHS)
      * \param a2  real array 2 (LHS)
      * \param b  shared real array to multiply both a1 and a2 (RHS)
      */
      void mulEqVPair(Array<double>& a1,
                      Array<double>& a2,
                      Array<double> const & b);

   } // namespace VecOp

} // namespace Pscf
#endif
