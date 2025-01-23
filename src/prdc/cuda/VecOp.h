#ifndef PRDC_VEC_OP_H
#define PRDC_VEC_OP_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "types.h"
#include <pscf/cuda/DeviceArray.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

/** 
* Functions that perform element-wise vector operations on the GPU.
*
* The functions defined in this file are wrappers for CUDA kernels that 
* perform the actual vector operations. The kernels themselves are only 
* intended to be called through their wrappers, so they are defined in 
* an anonymous namespace in VecOp.cu.
*
* The wrapper functions operate on DeviceArray objects containing 
* elements of type cudaReal or cudaComplex. The operations that are 
* performed by these functions include addition, subtraction, 
* multiplication, division, exponentiation, and assignment. The function 
* names will, correspondingly, begin with "add", "sub", "mul", "div", 
* "exp", or "eq" to indicate the operation being performed. Functions 
* are also included to perform compound assignment operations, i.e. 
* those that are performed using +=, -=, *=, and /= in C++. These 
* functions have names that begin with "addEq", "subEq", "mulEq", and 
* "divEq", respectively. 
*
* The functions are overloaded to perform their respective operations
* on any combination of cudaReal and cudaComplex input arrays, except
* those that would result in dividing by a complex number.
*
* The output (the LHS of the vector operation) will always be the first
* parameter passed to the function. The input argument(s) (on the RHS 
* of the vector operation) may be vectors or scalars. If an argument is 
* a vector (scalar), the function name will contain a V (S). For example, 
* addVV(A,B,C) implements vector-vector addition A[i] = B[i] + C[i], 
* while addVS(A,B,c) implements vector-scalar addition A[i] = B[i] + c 
* in which c is a scalar that is added to every element of B. In 
* operations involving both vectors and scalars, the vectors will always
* be listed first. So, for example, addVS exists, but addSV does not. 
* 
* Two wrapper functions are provided for each vector operation: 
* - The first accepts only the output array and the necessary input 
*   arrays / scalars. In these functions, each input array must be at 
*   least as long as the output array, and the element-wise operation 
*   will be performed for every element of the output array. All 
*   arrays will be indexed starting at element 0.
* - The second allows for vector operations to be performed using only
*   subsections (slices) of the input and output arrays. These functions 
*   require additional parameters: one index for each array involved in 
*   the operation (output and input), representing the element of each 
*   array at which to begin the slice, and an integer n, representing 
*   the size of the slices. Before calling the CUDA kernel, these 
*   functions check to ensure that the slices do not contain any indices
*   that exceed the length of the corresponding arrays.
* 
* Additional functions that perform multiple operations within a single 
* kernel are defined in VecOpMisc, which are not comprehensive and are
* written as-needed during the development of this software. VecOpMisc.h
* is included at the end of this file so that any code that includes
* this header file will also include VecOpMisc.h.
* 
* \ingroup Prdc_Cuda_Module 
* @{
*/
namespace VecOp {

// Assignment operations:
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void eqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
         const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void eqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  eqV(a, b, 0, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void eqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
         const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void eqV(DeviceArray<cudaComplex>& a, 
                DeviceArray<cudaComplex> const & b)
{  eqV(a, b, 0, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void eqS(DeviceArray<cudaReal>& a, cudaReal const b,
         const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void eqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  eqS(a, b, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void eqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
         const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void eqS(DeviceArray<cudaComplex>& a, cudaComplex const b)
{  eqS(a, b, 0, a.capacity()); }


// Addition operations
// ~~~~~~~~~~~~~~~~~~~

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void addVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaReal> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void addVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaComplex> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void addVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaComplex> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void addVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void addVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaReal> const & c)
{  addVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void addVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  cudaReal const c)
{  addVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaComplex const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void addVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaComplex const c)
{  addVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           cudaComplex const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
                  cudaComplex const c)
{  addVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void addVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaReal const c)
{  addVS(a, b, c, 0, 0, a.capacity()); }


// Subtraction operations
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void subVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaReal> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void subVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaComplex> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaComplex> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void subVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaReal> const & c)
{  subVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void subVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  cudaReal const c)
{  subVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaComplex const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void subVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaComplex const c)
{  subVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           cudaComplex const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void subVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaReal> const & b, 
                  cudaComplex const c)
{  subVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void subVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaReal const c)
{  subVS(a, b, c, 0, 0, a.capacity()); }


// Multiplication operations
// ~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void mulVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaReal> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void mulVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaComplex> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void mulVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaComplex> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void mulVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaReal> const & c)
{  mulVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void mulVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  cudaReal const c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaComplex const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void mulVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaComplex const c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           cudaComplex const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void mulVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaReal> const & b, cudaComplex const c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void mulVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaReal const c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }


// Division operations
// ~~~~~~~~~~~~~~~~~~~

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void divVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void divVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  DeviceArray<cudaReal> const & c)
{  divVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param beginIdC  index of the first entry to evaluate in array c
* \param n  the number of entries to evaluate
*/
void divVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaReal> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline void divVV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, 
                  DeviceArray<cudaReal> const & c)
{  divVV(a, b, c, 0, 0, 0, a.capacity()); }

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void divVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
                  cudaReal const c)
{  divVS(a, b, c, 0, 0, a.capacity()); }

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           cudaReal const c, const int beginIdA, 
           const int beginIdB, const int n);

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline void divVS(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b, cudaReal const c)
{  divVS(a, b, c, 0, 0, a.capacity()); }


// Exponentiation operations:
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void expV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
          const int beginIdA, const int beginIdB, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void expV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  expV(a, b, 0, 0, a.capacity()); }

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void expV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
          const int beginIdA, const int beginIdB, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void expV(DeviceArray<cudaComplex>& a, 
                 DeviceArray<cudaComplex> const & b)
{  expV(a, b, 0, 0, a.capacity()); }


// Compound operations: addition
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void addEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void addEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaComplex> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void addEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void addEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaReal> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void addEqS(DeviceArray<cudaReal>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void addEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void addEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
            const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void addEqS(DeviceArray<cudaComplex>& a, cudaComplex const b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void addEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void addEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
{  addEqS(a, b, 0, a.capacity()); }


// Compound operations: subtraction
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void subEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void subEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaComplex> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i]-=b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void subEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i]-=b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void subEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaReal> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void subEqS(DeviceArray<cudaReal>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void subEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void subEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
            const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void subEqS(DeviceArray<cudaComplex>& a, cudaComplex const b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void subEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void subEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
{  subEqS(a, b, 0, a.capacity()); }


// Compound operations: multiplication
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void mulEqV(DeviceArray<cudaReal>& a, 
                   DeviceArray<cudaReal> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void mulEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaComplex> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (mixed, b=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void mulEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (mixed, b=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void mulEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaReal> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void mulEqS(DeviceArray<cudaReal>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void mulEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void mulEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
            const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void mulEqS(DeviceArray<cudaComplex>& a, cudaComplex const b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i]*=b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void mulEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i]*=b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void mulEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
{  mulEqS(a, b, 0, a.capacity()); }


// Compound operations: division
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void divEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b)
{  divEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
void divEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline void divEqV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaReal> const & b)
{  divEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void divEqS(DeviceArray<cudaReal>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void divEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  divEqS(a, b, 0, a.capacity()); }

/**
* Vector division in-place, a[i] /= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
void divEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
            const int beginIdA, const int n);

/**
* Vector division in-place, a[i] /= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline void divEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
{  divEqS(a, b, 0, a.capacity()); }

/** @} */

} // namespace VecOp
} // namespace Cuda
} // namespace Prdc
} // namespace Pscf

#include "VecOpMisc.h" // Ensure that if VecOp is included, so is VecOpMisc

#endif
