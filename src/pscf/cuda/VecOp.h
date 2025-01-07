#ifndef PSCF_VEC_OP_H
#define PSCF_VEC_OP_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include "DeviceArray.h"
#include <util/containers/DArray.h>

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
* Additionally, a few miscellaneous functions are defined to perform 
* multiple of these operations within a single kernel, though these 
* are not overloaded to accommodate all input data types. In these 
* functions, we use "Vc" in the function name to denote that every
* element in the vector V is multiplied by a coefficient c. For 
* example, addEqVc represents the operation a[i] += b[i] * c. We also
* use "Many" to indicate that an undefined number of vectors (>2) are
* involved in an operation. For example, addVMany adds >2 vectors 
* together by passing an array of vectors, rather than a discrete set
* of vectors.
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
__host__ void eqV(DeviceArray<cudaReal>& a, 
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void eqV(DeviceArray<cudaReal>& a, 
                         DeviceArray<cudaReal> const & b)
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
__host__ void eqV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB, const int n);

/**
* Vector assignment, a[i] = b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void eqV(DeviceArray<cudaComplex>& a, 
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
__host__ void eqS(DeviceArray<cudaReal>& a, cudaReal const b,
                  const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void eqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  eqS(a, b, 0, a.capacity()); }

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void eqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                  const int beginIdA, const int n);

/**
* Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void eqS(DeviceArray<cudaComplex>& a, cudaComplex const b)
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
__host__ void addVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void addVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceArray<cudaComplex>& a, 
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
__host__ void addVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceArray<cudaComplex>& a, 
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
__host__ void addVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void addVV(DeviceArray<cudaComplex>& a, 
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
__host__ void addVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void addVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void addVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, cudaComplex const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void addVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void addVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void subVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void subVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceArray<cudaComplex>& a, 
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
__host__ void subVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void subVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void subVV(DeviceArray<cudaComplex>& a, 
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
__host__ void subVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void subVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void subVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, cudaComplex const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceArray<cudaComplex>& a, 
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
__host__ void subVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void subVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void mulVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void mulVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceArray<cudaComplex>& a, 
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
__host__ void mulVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceArray<cudaComplex>& a, 
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
__host__ void mulVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void mulVV(DeviceArray<cudaComplex>& a, 
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
__host__ void mulVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void mulVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void mulVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, cudaComplex const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void mulVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void mulVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void divVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void divVV(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void divVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
inline __host__ void divVV(DeviceArray<cudaComplex>& a, 
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
__host__ void divVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void divVS(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
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
__host__ void divVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n);

/**
* Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
inline __host__ void divVS(DeviceArray<cudaComplex>& a, 
                           DeviceArray<cudaComplex> const & b, 
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
__host__ void expV(DeviceArray<cudaReal>& a, 
                   DeviceArray<cudaReal> const & b,
                   const int beginIdA, const int beginIdB, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void expV(DeviceArray<cudaReal>& a, 
                          DeviceArray<cudaReal> const & b)
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
__host__ void expV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaComplex> const & b,
                   const int beginIdA, const int beginIdB, const int n);

/**
* Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void expV(DeviceArray<cudaComplex>& a, 
                          DeviceArray<cudaComplex> const & b)
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
__host__ void addEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void addEqV(DeviceArray<cudaReal>& a, 
                            DeviceArray<cudaReal> const & b)
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
__host__ void addEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void addEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void addEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector addition in-place, a[i] += b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void addEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void addEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void addEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  addEqS(a, b, 0, a.capacity()); }

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void addEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void addEqS(DeviceArray<cudaComplex>& a, 
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
__host__ void addEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector addition in-place, a[i] += b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void addEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
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
__host__ void subEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void subEqV(DeviceArray<cudaReal>& a, 
                            DeviceArray<cudaReal> const & b)
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
__host__ void subEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void subEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void subEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector subtraction in-place, a[i]-=b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void subEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void subEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void subEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  subEqS(a, b, 0, a.capacity()); }

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void subEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void subEqS(DeviceArray<cudaComplex>& a, 
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
__host__ void subEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void subEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
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
__host__ void mulEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void mulEqV(DeviceArray<cudaReal>& a, 
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
__host__ void mulEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void mulEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void mulEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (mixed, b=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void mulEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void mulEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void mulEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  mulEqS(a, b, 0, a.capacity()); }

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void mulEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void mulEqS(DeviceArray<cudaComplex>& a, 
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
__host__ void mulEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector multiplication in-place, a[i]*=b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void mulEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
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
__host__ void divEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void divEqV(DeviceArray<cudaReal>& a, 
                            DeviceArray<cudaReal> const & b)
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
__host__ void divEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n);

/**
* Vector division in-place, a[i] /= b[i], kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void divEqV(DeviceArray<cudaComplex>& a, 
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
__host__ void divEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void divEqS(DeviceArray<cudaReal>& a, cudaReal const b)
{  divEqS(a, b, 0, a.capacity()); }

/**
* Vector division in-place, a[i] /= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param beginIdA  index of the first entry to evaluate in array a
* \param n  the number of entries to evaluate
*/
__host__ void divEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n);

/**
* Vector division in-place, a[i] /= b, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
inline __host__ void divEqS(DeviceArray<cudaComplex>& a, cudaReal const b)
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
__host__ void addEqVc(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaReal> const & b, 
                      cudaReal const c, const int beginIdA, 
                      const int beginIdB, const int n);

/**
* Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
inline __host__ void addEqVc(DeviceArray<cudaReal>& a, 
                             DeviceArray<cudaReal> const & b, 
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
__host__ void sqNormV(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaComplex> const & b, 
                      const int beginIdA, const int beginIdB, const int n);

/**
* Squared norm of complex vector, a[i] = norm(b[i])^2, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
inline __host__ void sqNormV(DeviceArray<cudaReal>& a, 
                             DeviceArray<cudaComplex> const & b)
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
__host__ void expVc(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    cudaReal const c, const int beginIdA, 
                    const int beginIdB, const int n);

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
inline __host__ void expVc(DeviceArray<cudaReal>& a, 
                           DeviceArray<cudaReal> const & b, 
                           cudaReal const c)
{  expVc(a, b, c, 0, 0, a.capacity()); }

/**
* Add more than 2 vectors pointwise, GPU kernel.
*
* The input const pointer 'vecs' points to an array of const pointers. In
* other words, this is an array of arrays, where each array is represented
* by its pointer. The size of vecs is nVecs, and the size of vecs[i] is
* n (if i < nVecs). These nVecs vectors will be added and the result will
* be stored in vector 'a'.
*
* \param a  output array (LHS)
* \param vecs  array of pointers to DeviceArrays to be added
* \param nVecs  number of vectors to be added
* \param n  size of arrays
*/
__global__ void _addVMany(cudaReal* a, cudaReal const ** vecs,
                          const int nVecs, const int n);

/**
* Add more than 2 vectors pointwise, kernel wrapper.
*
* The input array 'vecs' contains const pointers to each array that will be 
* added together. The size of vecs determines the number of vectors that
* will ultimately be added together by the GPU kernel.
*
* \param a  output array (LHS)
* \param vecs  array of pointers to DeviceArrays to be added
*/
__host__ void addVMany(DeviceArray<cudaReal>& a, 
                       DArray<DeviceArray<cudaReal> const *> const & vecs);

/**
* Multiply more than 2 vectors pointwise, GPU kernel.
*
* The input const pointer 'vecs' points to an array of const pointers. In
* other words, this is an array of arrays, where each array is represented
* by its pointer. The size of vecs is nVecs, and the size of vecs[i] is
* n (if i < nVecs). These nVecs vectors will be multiplied and the result 
* will be stored in vector 'a'.
*
* \param a  output array (LHS)
* \param vecs  array of pointers to DeviceArrays to be multiplied
* \param nVecs  number of vectors to be multiplied
* \param n  size of arrays
*/
__global__ void _mulVMany(cudaReal* a, cudaReal const ** vecs,
                          const int nVecs, const int n);

/**
* Multiply more than 2 vectors pointwise, kernel wrapper.
*
* The input array 'vecs' contains const pointers to each array that will be 
* multiplied together. The size of vecs determines the number of vectors
* that will ultimately be multiplied together by the GPU kernel.
*
* \param a  output array (LHS)
* \param vecs  array of pointers to DeviceArrays to be multiplied
*/
__host__ void mulVMany(DeviceArray<cudaReal>& a, 
                       DArray<DeviceArray<cudaReal> const *> const & vecs);

/** @} */

} // namespace VecOp
} // namespace Pscf
#endif
