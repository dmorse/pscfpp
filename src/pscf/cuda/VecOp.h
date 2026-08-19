#ifndef PSCF_CUDA_VEC_OP_H
#define PSCF_CUDA_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>

namespace Pscf {

   /**
   * Functions that perform element-wise vector operations on the GPU.
   *
   * The functions declared in this header operate on DeviceArray objects
   * with elements of type cudaReal or cudaComplex. The operations that
   * are performed by these functions include assignment, addition,
   * subtraction, multiplication, division, exponentiation, square and
   * absolute magnitude. Function names, correspondingly, begin with 
   * "eq", "add", "sub", "mul", "div", "exp", "sq" and "abs" to indicate
   * the operation being performed. Functions that perform in-place 
   * arithmetic assignment operations, which are analogous to those 
   * performed using the +=, -=, *=, and /= C/C++ operators, have names 
   * that begin with "addEq", "subEq", "mulEq", and "divEq".
   *
   * These functions are overloaded to perform their respective operations
   * on any combination of cudaReal and cudaComplex input arrays, except
   * those that would result in division by a complex number.
   *
   * The output (the LHS of the vector operation) is always the first
   * parameter passed to the function. The input argument(s) (on the RHS
   * of the vector operation) may be vectors or scalars. If an argument
   * is a vector (scalar), the function name will contain a V (S). For
   * example, addVV(A,B,C) implements the vector-vector addition operation
   * A[i] = B[i] + C[i], while addVS(A,B,c) implements vector-scalar
   * addition A[i] = B[i] + c, in which c is a scalar that is added to
   * every element of array B. In commutative binary operations involving
   * a vector and a scalar, the vector is listed first. So, for example,
   * addVS exists, but addSV does not.
   *
   * Two wrapper functions are provided for each vector operation:
   * - One function accepts only an output array and any required input
   *   arrays and scalars. In this function, each input array must be at
   *   least as long as the output array, and the element-wise operation
   *   is performed for every element of the output array, starting from
   *   element 0 of all input and output arrays.
   * - One function allows for vector operations to be performed using
   *   subsections (slices) of the input and output arrays. This function
   *   require additional parameters: one index for each input and output
   *   array to indicate the index of the first element of the slice in
   *   that array, and an integer n indicating the size of the slice for
   *   all arrays. Before calling the CUDA kernel, such functions check
   *   to ensure that none of the slices exceeds the array bounds.
   *
   * Additional functions that perform multiple operations within a
   * single kernel are defined in VecOpMisc. This collection is not
   * comprehensive and is added to as-needed during the development of
   * this software.  VecOpMisc.h is included at the end of VecOp.h so
   * that any code that includes VecOp.h will also include VecOpMisc.h.
   *
   * The functions declared in this file are wrappers for CUDA kernels
   * that perform the actual parallel vector operations. The underlying
   * CUDA kernels are only intended to be called through their wrappers.
   * To enforce this, kernels are defined in an anonymous namespace in 
   * the file VecOp.cu, and are thus only accessible for use within that 
   * source file.
   *
   * \ingroup Pscf_Cuda_Module
   * \defgroup Pscf_Cuda_VecOp_Module VecOp (GPU)
   */

   /**
   *  Vector operations on GPU or CPU.
   * 
   * \ingroup Pscf_Math_Module
   */
   namespace VecOp {

      // Assignment operations:

      /**
      * Vector-vector assignment, a[i] = b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void eqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector assignment, a[i] = b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void eqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b)
      {  eqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector assignment, a[i] = b[i], (real, device to host).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void eqV(Array<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector assignment, a[i] = b[i], (real, device to host).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void eqV(Array<cudaReal>& a,
               DeviceArray<cudaReal> const & b)
      {  eqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector assignment, a[i] = b[i], (real, host to device).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void eqV(DeviceArray<cudaReal>& a,
               Array<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector assignment, a[i] = b[i], (real, host to device).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void eqV(DeviceArray<cudaReal>& a,
               Array<cudaReal> const & b)
      {  eqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector assignment, a[i] = b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void eqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector assignment, a[i] = b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      inline
      void eqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b)
      {  eqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-scalar assignment, a[i] = b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void eqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n);

      /**
      * Vector-scalar assignment, a[i] = b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void eqS(DeviceArray<cudaReal>& a, const cudaReal b)
      {  eqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar assignment, a[i] = b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void eqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n);

      /**
      * Vector-scalar assignment, a[i] = b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      inline
      void eqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
      {  eqS(a, b, 0, a.capacity()); }

      // Addition operations

      /**
      * Vector-vector addition, a[i] = b[i] + c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void addVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector addition, a[i] = b[i] + c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void addVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c)
      {  addVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector addition, a[i] = b[i] + c[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void addVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaComplex> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector addition, a[i] = b[i] + c[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      inline
      void addVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaComplex> const & c)
      {  addVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void addVV(DeviceArray<cudaComplex> & a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaComplex> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex array (RHS)
      */
      inline
      void addVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaComplex> const & c)
      {  addVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void addVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector addition, a[i] = b[i] + c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void addVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c)
      {  addVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-scalar addition, a[i] = b[i] + c (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB,
                 const int n);

      /**
      * Vector-scalar addition, a[i] = b[i] + c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void addVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c)
      {  addVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar addition, a[i] = b[i] + c, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaComplex c,
                 const int beginIdA, const int beginIdB,
                 const int n);

      /**
      * Vector addition, a[i] = b[i] + c, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      inline
      void addVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaComplex c)
      {  addVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector addition, a[i] = b[i] + c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaComplex c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector addition, a[i] = b[i] + c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      */
      inline
      void addVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaComplex c)
      {  addVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void addVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c)
      {  addVS(a, b, c, 0, 0, a.capacity()); }


      // Subtraction operations

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void subVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c,
      	   const int beginIdA, const int beginIdB, const int beginIdC,
      	   const int n);

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void subVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c)
      {  subVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void subVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaComplex> const & c, const int beginIdA,
                 const int beginIdB, const int beginIdC, const int n);

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      inline
      void subVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaComplex> const & c)
      {  subVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void subVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaComplex> const & c,
                 const int beginIdA,
                 const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex array (RHS)
      */
      inline
      void subVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaComplex> const & c)
      {  subVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (mixed, c = real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void subVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA,
                 const int beginIdB, const int beginIdC, const int n);

      /**
      * Vector-vector subtraction, a[i] = b[i] - c[i] (mixed, c = real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real  array (RHS)
      */
      inline
      void subVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c)
      {  subVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB,
                 const int n);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void subVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c)
      {  subVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaComplex c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      inline
      void subVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaComplex c)
      {  subVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector subtraction, a[i] = b[i] - c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaComplex c,
      	   const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector subtraction, a[i] = b[i] - c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      */
      inline
      void subVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaComplex c)
      {  subVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector subtraction, a[i] = b[i] - c (mixed, c = real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector subtraction, a[i] = b[i] - c (mixed, c = real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void subVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c)
      {  subVS(a, b, c, 0, 0, a.capacity()); }


      // Multiplication operations

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void mulVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c,
      	   const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void mulVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c)
      {  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void mulVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaComplex> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      inline
      void mulVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaComplex> const & c)
      {  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector multiplication, a[i] = b[i] * c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void mulVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaComplex> const & c, const int beginIdA,
                 const int beginIdB, const int beginIdC, const int n);

      /**
      * Vector-vector multiplication, a[i]=b[i]*c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      */
      inline
      void mulVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaComplex> const & c)
      {  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector multiplication, a[i]=b[i]*c[i] (mixed, c = real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void mulVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector multiplication, a[i]=b[i]*c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void mulVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c)
      {  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB,
                 const int n);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void mulVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c)
      {  mulVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaComplex c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      inline
      void mulVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
      	   const cudaComplex c)
      {  mulVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaComplex c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      */
      inline
      void mulVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaReal> const & b,
      	   const cudaComplex c)
      {  mulVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void mulVS(DeviceArray<cudaComplex>& a,
                        DeviceArray<cudaComplex> const & b,
                        const cudaReal c)
      {  mulVS(a, b, c, 0, 0, a.capacity()); }


      // Division operations

      /**
      * Vector-vector division, a[i] = b[i] / c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void divVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector division, a[i] = b[i] / c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void divVV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 DeviceArray<cudaReal> const & c)
      {  divVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void divVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA, const int beginIdB, const int beginIdC,
                 const int n);

      /**
      * Vector-vector division, a[i] = b[i] / c[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      inline
      void divVV(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 DeviceArray<cudaReal> const & c)
      {  divVV(a, b, c, 0, 0, 0, a.capacity()); }

      /**
      * Vector-scalar division, a[i] = b[i] / c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void divVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c, const int beginIdA,
                 const int beginIdB, const int n);

      /**
      * Vector-scalar division, a[i] = b[i] / c, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void divVS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b,
                 const cudaReal c)
      {  divVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Vector-scalar division, a[i] = b[i] / c (mixed, c real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void divVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c,
                 const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-scalar division, a[i] = b[i] / c (mixed, c = real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      inline
      void divVS(DeviceArray<cudaComplex>& a,
                 DeviceArray<cudaComplex> const & b,
                 const cudaReal c)
      {  divVS(a, b, c, 0, 0, a.capacity()); }

      /**
      * Scalar-vector division, a[i] = b / c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param c  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdC  index of first element in a slice of array c
      * \param n  number of elements in the slice
      */
      void divSV(DeviceArray<cudaReal>& a,
                 const cudaReal b,
                 DeviceArray<cudaReal> const & c,
                 const int beginIdA, const int beginIdC, const int n);

      /**
      * Scalar-vector division, a[i] = b / c[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param c  real array (RHS)
      */
      inline
      void divSV(DeviceArray<cudaReal>& a,
                 const cudaReal b,
                 DeviceArray<cudaReal> const & c)
      {  divSV(a, b, c, 0, 0, a.capacity()); }

      // In-place addition

      /**
      * Vector in-place addition, a[i] += b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place addition, a[i] += b[i] (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void addEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b)
      {  addEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector in-place addition, a[i] += b[i] (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place addition, a[i] += b[i] (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      inline
      void addEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaComplex> const & b)
      {  addEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector in-place addition, a[i] += b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void addEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place addition, a[i] += b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void addEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b)
      {  addEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-scalar in-place addition, a[i] += b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void addEqS(DeviceArray<cudaReal>& a,
                  const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place addition, a[i] += b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void addEqS(DeviceArray<cudaReal>& a, const cudaReal b)
      {  addEqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar in-place addition, a[i] += b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void addEqS(DeviceArray<cudaComplex>& a,
                  const cudaComplex b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place addition, a[i] += b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      */
      inline
      void addEqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
      {  addEqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar in-place addition, a[i] += b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void addEqS(DeviceArray<cudaComplex>& a,
                  const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place addition, a[i] += b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void addEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
      {  addEqS(a, b, 0, a.capacity()); }


      // In-place subtraction

      /**
      * Vector in-place subtraction, a[i] -= b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector in-place subtraction, a[i] -= b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void subEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b)
      {  subEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector in-place subtraction, a[i] -= b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector in-place subtraction, a[i] -= b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      inline
      void subEqV(DeviceArray<cudaComplex>& a,
                         DeviceArray<cudaComplex> const & b)
      {  subEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector in-place subtraction, a[i]-=b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void subEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector in-place subtraction, a[i]-=b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      */
      inline
      void subEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b)
      {  subEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector in-place subtraction, a[i] -= b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void subEqS(DeviceArray<cudaReal>& a, const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector in-place subtraction, a[i] -= b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void subEqS(DeviceArray<cudaReal>& a, const cudaReal b)
      {  subEqS(a, b, 0, a.capacity()); }

      /**
      * Vector in-place subtraction, a[i] -= b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void subEqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
                  const int beginIdA, const int n);

      /**
      * Vector in-place subtraction, a[i] -= b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      inline
      void subEqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
      {  subEqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void subEqS(DeviceArray<cudaComplex>& a,
                  const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void subEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
      {  subEqS(a, b, 0, a.capacity()); }


      // In-place multiplication

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void mulEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b)
      {  mulEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      inline
      void mulEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaComplex> const & b)
      {  mulEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector in-place multiplication, a[i]*=b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void mulEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place multiplication, a[i] *= b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void mulEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b)
      {  mulEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-scalar in-place multiplication, a[i] *= b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void mulEqS(DeviceArray<cudaReal>& a, const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void mulEqS(DeviceArray<cudaReal>& a, const cudaReal b)
      {  mulEqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar in-place multiplication, a[i] *= b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void mulEqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b, (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      inline
      void mulEqS(DeviceArray<cudaComplex>& a, const cudaComplex b)
      {  mulEqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar in-place multiplication, a[i]*=b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void mulEqS(DeviceArray<cudaComplex>& a,
                  const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place multiplication, a[i]*=b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void mulEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
      {  mulEqS(a, b, 0, a.capacity()); }


      // In-place division

      /**
      * Vector-vector in-place division, a[i] /= b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void divEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place division, a[i] /= b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void divEqV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b)
      {  divEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-vector in-place division, a[i] /= b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void divEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

      /**
      * Vector-vector in-place division, a[i] /= b[i] (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void divEqV(DeviceArray<cudaComplex>& a,
                  DeviceArray<cudaReal> const & b)
      {  divEqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector-scalar in-place division, a[i] /= b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void divEqS(DeviceArray<cudaReal>& a,
                  const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place division, a[i] /= b, (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void divEqS(DeviceArray<cudaReal>& a, const cudaReal b)
      {  divEqS(a, b, 0, a.capacity()); }

      /**
      * Vector-scalar in-place division, a[i] /= b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param n  number of elements in the slice
      */
      void divEqS(DeviceArray<cudaComplex>& a,
                  const cudaReal b,
                  const int beginIdA, const int n);

      /**
      * Vector-scalar in-place division, a[i] /= b (mixed).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      inline
      void divEqS(DeviceArray<cudaComplex>& a, const cudaReal b)
      {  divEqS(a, b, 0, a.capacity()); }

      // Exponentiation operations

      /**
      * Vector exponentiation, a[i] = exp(b[i]), (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void expV(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b,
                const int beginIdA, const int beginIdB,
                const int n);

      /**
      * Vector exponentiation, a[i] = exp(b[i]), (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void expV(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b)
      {  expV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector exponentiation, a[i] = exp(b[i]), (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void expV(DeviceArray<cudaComplex>& a,
                DeviceArray<cudaComplex> const & b,
                const int beginIdA, const int beginIdB,
                const int n);

      /**
      * Vector exponentiation, a[i] = exp(b[i]), (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      inline
      void expV(DeviceArray<cudaComplex>& a,
                DeviceArray<cudaComplex> const & b)
      {  expV(a, b, 0, 0, a.capacity()); }

      // Vector (element-wise) square

      /**
      * Vector square, a[i] = b[i]*b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void sqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB,
               const int n);

      /**
      * Vector square, a[i] = b[i]*b[i], (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void sqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b)
      {  sqV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector square, a[i] = b[i]*b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  number of elements in the slice
      */
      void sqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB,
               const int n);

      /**
      * Vector square, a[i] = b[i]*b[i], (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      inline
      void sqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b)
      {  sqV(a, b, 0, 0, a.capacity()); }

      // Absolute magnitude

      /**
      * Vector absolute magnitude, a[i] = abs(b[i]) (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  size of arrays
      */
      void absV(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b,
                const int beginIdA, const int beginIdB,
                const int n);

      /**
      * Vector absolute magnitude, a[i] = abs(b[i]) (real).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      */
      inline
      void absV(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b)
      {  absV(a, b, 0, 0, a.capacity()); }

      /**
      * Vector squared absolute magnitude, a[i] = |b[i]|^2 (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      * \param beginIdA  index of first element in a slice of array a
      * \param beginIdB  index of first element in a slice of array b
      * \param n  size of arrays
      */
      void sqAbsV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB,
                  const int n);

      /**
      * Vector absolute magnitude squared, a[i] = |b[i]|^2 (complex).
      *
      * \ingroup Pscf_Cuda_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  conplex array (RHS)
      */
      inline
      void sqAbsV(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaComplex> const & b)
      {  sqAbsV(a, b, 0, 0, a.capacity()); }

   } // namespace VecOp
} // namespace Pscf

// Ensure that if VecOp.h is included, so is VecOpMisc.h
#include "VecOpMisc.h"
#endif
