#ifndef PSCF_VEC_OP_H
#define PSCF_VEC_OP_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include "DeviceDArray.h"

namespace Pscf {
namespace VecOp {

/** 
* Functions that perform element-wise vector operations on the GPU.
*
* CUDA kernels that perform the vector operations are defined here, as 
* well as wrapper functions to be called by the host CPU, which call 
* the kernel internally. The kernels should always be accessed through 
* the wrapper functions when possible. 
*
* The wrapper functions operate on DeviceDArray objects containing 
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
* Additionally, a few miscellaneous functions are defined to perform 
* multiple of these operations within a single kernel, though these 
* are not overloaded to accommodate all input data types. In these 
* functions, we use "Vc" in the function name to denote that every
* element in the vector V is multiplied by a coefficient c. For 
* example, addEqVc represents the operation a[i] += b[i] * c.
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
* \ingroup Pscf_Cuda_Module 
* @{
*/


// Assignment operations:
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector assignment, a[i] = b[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _eqV(cudaReal* a, cudaReal const * b, const int n);

/**
* Vector assignment, a[i] = b[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _eqV(cudaComplex* a, cudaComplex const * b, const int n);

/**
* Vector assignment, a[i] = b, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _eqS(cudaReal* a, cudaReal const b, const int n);

/**
* Vector assignment, a[i] = b, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _eqS(cudaComplex* a, cudaComplex const b, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void eqV(DeviceDArray<cudaReal>& a, 
                  DeviceDArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void eqV(DeviceDArray<cudaReal>& a, 
                         DeviceDArray<cudaReal> const & b)
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
__host__ void eqV(DeviceDArray<cudaComplex>& a, 
                  DeviceDArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void eqV(DeviceDArray<cudaComplex>& a, 
                         DeviceDArray<cudaComplex> const & b)
{  eqV(a, b, 0, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void eqS(DeviceDArray<cudaReal>& a, cudaReal const b,
                  const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void eqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{  eqS(a, b, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void eqS(DeviceDArray<cudaComplex>& a, cudaComplex const b,
                  const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void eqS(DeviceDArray<cudaComplex>& a, cudaComplex const b)
{  eqS(a, b, 0, a.capacity()); }


// Addition operations
// ~~~~~~~~~~~~~~~~~~~

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addVV(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const c, const int n);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

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
__host__ void addVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaComplex> const & c)
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
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaComplex> const & c)
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
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void addVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
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
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaComplex const c)
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
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, cudaComplex const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaReal> const & b, 
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
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaReal const c)
{  addVS(a, b, c, 0, 0, a.capacity()); }


// Subtraction operations
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subVV(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

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
__host__ void subVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaComplex> const & c)
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
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaComplex> const & c)
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
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void subVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
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
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaComplex const c)
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
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, cudaComplex const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaReal> const & b, 
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
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaReal const c)
{  subVS(a, b, c, 0, 0, a.capacity()); }


// Multiplication operations
// ~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector multiplication, a[i] = b[i] * c[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulVV(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulVS(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

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
__host__ void mulVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaComplex> const & c)
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
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaComplex> const & c)
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
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void mulVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
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
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaComplex const c)
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
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, cudaComplex const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           cudaComplex const c)
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
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaReal const c)
{  mulVS(a, b, c, 0, 0, a.capacity()); }


// Division operations
// ~~~~~~~~~~~~~~~~~~~

/**
* Vector division, a[i] = b[i] / c[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _divVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n);

/**
* Vector division, a[i] = b[i] / c[i], GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _divVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n);

/**
* Vector division, a[i] = b[i] / c, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n);

/**
* Vector division, a[i] = b[i] / c, GPU kernel (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

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
__host__ void divVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void divVV(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void divVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void divVV(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           DeviceDArray<cudaReal> const & c)
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
__host__ void divVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void divVS(DeviceDArray<cudaReal>& a, 
                           DeviceDArray<cudaReal> const & b, 
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
__host__ void divVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void divVS(DeviceDArray<cudaComplex>& a, 
                           DeviceDArray<cudaComplex> const & b, 
                           cudaReal const c)
{  divVS(a, b, c, 0, 0, a.capacity()); }


// Exponentiation operations:
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector exponentiation, a[i] = exp(b[i]), GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _expV(cudaReal* a, cudaReal const * b, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _expV(cudaComplex* a, cudaComplex const * b, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void expV(DeviceDArray<cudaReal>& a, 
                   DeviceDArray<cudaReal> const & b,
                   const int beginIdA, const int beginIdB, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void expV(DeviceDArray<cudaReal>& a, 
                          DeviceDArray<cudaReal> const & b)
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
__host__ void expV(DeviceDArray<cudaComplex>& a, 
                   DeviceDArray<cudaComplex> const & b,
                   const int beginIdA, const int beginIdB, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void expV(DeviceDArray<cudaComplex>& a, 
                          DeviceDArray<cudaComplex> const & b)
{  expV(a, b, 0, 0, a.capacity()); }


// Compound operations: addition
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector addition in-place, a[i] += b[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addEqV(cudaReal* a, cudaReal const * b, const int n);

/**
* Vector addition in-place, a[i] += b[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addEqV(cudaComplex* a, cudaComplex const * b, const int n);

/**
* Vector addition in-place, a[i] += b[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addEqV(cudaComplex* a, cudaReal const * b, const int n);

/**
* Vector addition in-place, a[i] += b, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addEqS(cudaReal* a, cudaReal const b, const int n);

/**
* Vector addition in-place, a[i] += b, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addEqS(cudaComplex* a, cudaComplex const b, const int n);

/**
* Vector addition in-place, a[i] += b, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addEqS(cudaComplex* a, cudaReal const b, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void addEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void addEqV(DeviceDArray<cudaReal>& a, 
                            DeviceDArray<cudaReal> const & b)
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
__host__ void addEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void addEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaComplex> const & b)
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
__host__ void addEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void addEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaReal> const & b)
{  addEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void addEqS(DeviceDArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void addEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void addEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void addEqS(DeviceDArray<cudaComplex>& a, 
                            cudaComplex const b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void addEqS(DeviceDArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void addEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{  addEqS(a, b, 0, a.capacity()); }


// Compound operations: subtraction
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector subtraction in-place, a[i] -= b[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subEqV(cudaReal* a, cudaReal const * b, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subEqV(cudaComplex* a, cudaComplex const * b, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subEqV(cudaComplex* a, cudaReal const * b, const int n);

/**
* Vector subtraction in-place, a[i] -= b, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subEqS(cudaReal* a, cudaReal const b, const int n);

/**
* Vector subtraction in-place, a[i] -= b, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subEqS(cudaComplex* a, cudaComplex const b, const int n);

/**
* Vector subtraction in-place, a[i] -= b, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subEqS(cudaComplex* a, cudaReal const b, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void subEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void subEqV(DeviceDArray<cudaReal>& a, 
                            DeviceDArray<cudaReal> const & b)
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
__host__ void subEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void subEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaComplex> const & b)
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
__host__ void subEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i]-=b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void subEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaReal> const & b)
{  subEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void subEqS(DeviceDArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void subEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void subEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void subEqS(DeviceDArray<cudaComplex>& a, 
                            cudaComplex const b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void subEqS(DeviceDArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void subEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{  subEqS(a, b, 0, a.capacity()); }


// Compound operations: multiplication
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector multiplication in-place, a[i] *= b[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulEqV(cudaReal* a, cudaReal const * b, const int n);

/**
* Vector multiplication in-place, a[i] *= b[i], GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulEqV(cudaComplex* a, cudaComplex const * b, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulEqV(cudaComplex* a, cudaReal const * b, const int n);

/**
* Vector multiplication in-place, a[i] *= b, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulEqS(cudaReal* a, cudaReal const b, const int n);

/**
* Vector multiplication in-place, a[i] *= b, GPU kernel (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulEqS(cudaComplex* a, cudaComplex const b, const int n);

/**
* Vector multiplication in-place, a[i] *= b, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulEqS(cudaComplex* a, cudaReal const b, const int n);

/**
* Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void mulEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void mulEqV(DeviceDArray<cudaReal>& a, 
                            DeviceDArray<cudaReal> const & b)
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
__host__ void mulEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void mulEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaComplex> const & b)
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
__host__ void mulEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (mixed, b=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void mulEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaReal> const & b)
{  mulEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void mulEqS(DeviceDArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void mulEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void mulEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void mulEqS(DeviceDArray<cudaComplex>& a, 
                            cudaComplex const b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i]*=b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void mulEqS(DeviceDArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i]*=b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void mulEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{  mulEqS(a, b, 0, a.capacity()); }


// Compound operations: division
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector division in-place, a[i] /= b[i], GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _divEqV(cudaReal* a, cudaReal const * b, const int n);

/**
* Vector division in-place, a[i] /= b[i], GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _divEqV(cudaComplex* a, cudaReal const * b, const int n);

/**
* Vector division in-place, a[i] /= b, GPU kernel (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divEqS(cudaReal* a, cudaReal const b, const int n);

/**
* Vector division in-place, a[i] /= b, GPU kernel (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divEqS(cudaComplex* a, cudaReal const b, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void divEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void divEqV(DeviceDArray<cudaReal>& a, 
                            DeviceDArray<cudaReal> const & b)
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
__host__ void divEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void divEqV(DeviceDArray<cudaComplex>& a, 
                            DeviceDArray<cudaReal> const & b)
{  divEqV(a, b, 0, 0, a.capacity()); }

/**
* Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void divEqS(DeviceDArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void divEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{  divEqS(a, b, 0, a.capacity()); }

/**
* Vector division in-place, a[i] /= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void divEqS(DeviceDArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector division in-place, a[i] /= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void divEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{  divEqS(a, b, 0, a.capacity()); }


// Miscellaneous combined operations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector addition in-place w/ coefficient, a[i] += b[i] * c, GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
* \param n  size of arrays
*/
__global__ void _addEqVc(cudaReal* a, cudaReal const * b, 
                            cudaReal const c, const int n);

/**
* Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void addEqVc(DeviceDArray<cudaReal>& a, 
                         DeviceDArray<cudaReal> const & b, 
                         cudaReal const c, const int beginIdA, 
                         const int beginIdB, const int n);

/**
* Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
inline __host__ void addEqVc(DeviceDArray<cudaReal>& a, 
                                DeviceDArray<cudaReal> const & b, 
                                cudaReal const c)
{  addEqVc(a, b, c, 0, 0, a.capacity()); }

/**
* Squared norm of complex vector, a[i] = norm(b[i])^2, GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _sqNormV(cudaReal* a, cudaComplex const * b, const int n);

/**
* Squared norm of complex vector, a[i] = norm(b[i])^2, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void sqNormV(DeviceDArray<cudaReal>& a, 
                      DeviceDArray<cudaComplex> const & b, 
                      const int beginIdA, const int beginIdB, const int n);

/**
* Squared norm of complex vector, a[i] = norm(b[i])^2, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void sqNormV(DeviceDArray<cudaReal>& a, 
                             DeviceDArray<cudaComplex> const & b)
{  sqNormV(a, b, 0, 0, a.capacity()); }

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
* \param n  size of arrays
*/
__global__ void _expVc(cudaReal* a, cudaReal const * b, 
                          cudaReal const c, const int n);

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
* \param beginIdA  index of the first entry to evaluate in array a
* \param beginIdB  index of the first entry to evaluate in array b
* \param n  the number of entries to evaluate
*/
__host__ void expVc(DeviceDArray<cudaReal>& a, 
                       DeviceDArray<cudaReal> const & b, 
                       cudaReal const c, const int beginIdA, 
                       const int beginIdB, const int n);

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
inline __host__ void expVc(DeviceDArray<cudaReal>& a, 
                              DeviceDArray<cudaReal> const & b, 
                              cudaReal const c)
{  expVc(a, b, c, 0, 0, a.capacity()); }

/** @} */

} // namespace VecOp
} // namespace Pscf
#endif
