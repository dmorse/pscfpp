#ifndef RPC_VEC_OP_H
#define RPC_VEC_OP_H

#include <util/containers/Array.h>

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /** 
   * Functions that perform element-wise vector operations on the Cpu.
   *
   * Operations that are performed by these functions include addition, 
   * subtraction, multiplication, division, exponentiation, and assignment. 
   * The function names will, correspondingly, begin with "add", "sub", 
   * "mul", "div", "exp", or "eq" to indicate the relevant operation.
   * Functions are also included to perform compound assignment operations, 
   * i.e.  those that are performed using +=, -=, *=, and /= in C++. These 
   * functions have names that begin with "addEq", "subEq", "mulEq", and 
   * "divEq", respectively. 
   *
   * The output (the LHS of the vector operation) is always the first
   * parameter passed to the function. The input argument(s) (on the RHS 
   * of the vector operation) may be vectors or scalars. If an argument is 
   * a vector (scalar), the function name will contain a V (S). For 
   * example, the function addVV(A,B,C) implements vector-vector addition 
   * A[i] = B[i] + C[i], while addVS(A,B,c) implements vector-scalar 
   * addition A[i] = B[i] + c in which c is a scalar that is added to every
   *  element of B. In commutative operations involving both vectors and 
   * scalars, the vectors are listed first. So, for example, addVS exists, 
   * but addSV does not. 
   * 
   * \ingroup Prdc_Cpu_Module 
   * @{
   */
   namespace VecOp {
   
      // Assignment
      
      /**
      * Vector assignment, a[i] = b[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void eqV(Array<double>& a, Array<double> const & b);
      
      /**
      * Vector assignment, a[i] = b.
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void eqS(Array<double>& a, double b);
      
      // Addition
      
      /**
      * Vector addition, a[i] = b[i] + c[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void addVV(Array<double>& a, Array<double> const & b, 
                 Array<double> const & c);
      
      /**
      * Vector addition, a[i] = b[i] + c.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void addVS(Array<double>& a, Array<double> const & b, double c);
      
      // Subtraction
      
      /**
      * Vector subtraction, a[i] = b[i] - c[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void subVV(Array<double>& a, 
                 Array<double> const & b, Array<double> const & c);
      
      /**
      * Vector subtraction, a[i] = b[i] - c.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void subVS(Array<double>& a, Array<double> const & b, double c);
      
      // Multiplication
      
      /**
      * Vector multiplication, a[i] = b[i] * c[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void mulVV(Array<double>& a, 
                 Array<double> const & b, Array<double> const & c);
      
      /**
      * Vector multiplication, a[i] = b[i] * c.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void mulVS(Array<double>& a, Array<double> const & b, double c);
      
      // Division
      
      /**
      * Vector division, a[i] = b[i] / c[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      void divVV(Array<double>& a, 
                 Array<double> const & b, Array<double> const & c);
      
      /**
      * Vector division, a[i] = b[i] / c.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      */
      void divVS(Array<double>& a, Array<double> const & b, double c);
      
      /**
      * Vector division, a[i] = b / c[i].
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param c  input array (RHS)
      */
      void divSV(Array<double>& a, double b, Array<double> const & c);
      
      // Exponentiation
      
      /**
      * Vector exponentiation, a[i] = exp(b[i]).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void expV(Array<double>& a, Array<double> const & b);
      
      
      // Compound addition
      
      /**
      * Vector addition in-place, a[i] += b[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void addEqV(Array<double>& a, Array<double> const & b);
      
      /**
      * Vector addition in-place, a[i] += b.
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void addEqS(Array<double>& a, double b);
      
      // Compound subtraction
      
      /**
      * Vector subtraction in-place, a[i] -= b[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void subEqV(Array<double>& a, Array<double> const & b);
      
      /**
      * Vector subtraction in-place, a[i] -= b.
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void subEqS(Array<double>& a, double b);
      
      // Compound multiplication
      
      /**
      * Vector multiplication in-place, a[i] *= b[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void mulEqV(Array<double>& a, Array<double> const & b);
      
      /**
      * Vector multiplication in-place, a[i] *= b.
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void mulEqS(Array<double>& a, double b);
      
      // Compound division
      
      /**
      * Vector division in-place, a[i] /= b[i].
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      void divEqV(Array<double>& a, Array<double> const & b);
      
      /**
      * Vector division in-place, a[i] /= b.
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      void divEqS(Array<double>& a, double b);
   
   /** @} */

   } // namespace VecOp

} // namespace Cpu
} // namespace Prdc
} // namespace Pscf

#endif
