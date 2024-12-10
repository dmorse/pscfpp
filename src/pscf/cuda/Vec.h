#ifndef PSCF_VEC_H
#define PSCF_VEC_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include "DeviceDArray.h"

namespace Pscf {
namespace Vec {

/** 
* Functions that perform element-wise vector operations on the GPU.
*
* The operations that are performed by these functions include addition,
* subtraction, multiplication, division, exponentiation, and assignment. 
* The function names will, correspondingly, begin with "add", "sub", 
* "mul", "div", "exp", or "eq" to indicate the operation being performed. 
* Functions are also included to perform compound assignment operations, 
* i.e. those that are performed using +=, -=, *=, and /= in C++. These 
* functions have names that begin with "addEq", "subEq", "mulEq", and 
* "divEq", respectively. Finally, some functions are defined to perform
* multiple of these operations within a single kernel, though these are
* not comprehensive and are defined only as-needed within the software. 
*
* Each function will be overloaded to accept multiple data types. The
* output (the LHS of the vector operation) will be a vector of type 
* cudaReal or cudaComplex, and will always be the first parameter 
* passed to the function. The input argument(s) (on the RHS of the 
* vector operation) may be vectors or scalars of type cudaReal or 
* cudaComplex. Any combination of inputs that leads to an output vector 
* of type cudaReal or cudaComplex are included among the set of 
* overloaded functions. All operations are performed element-wise.
*
* If an argument is a vector (scalar), the function name will contain 
* a V (S). For example, addVV(A,B,C) implements vector-vector addition 
* A[i] = B[i] + C[i], while addVS(A,B,c) implements vector-scalar 
* addition A[i] = B[i] + c in which c is a scalar that is added to 
* every element of B.
*
* In operations involving both vectors and scalars, the vectors will
* always be listed first. So, for example, addVS exists, but addSV does 
* not. In operations involving two vectors of different types where the
* order does not matter (addition and multiplication), the cudaReal 
* vector will always be listed first. We do not include any functions
* that divide by a complex number.
* 
* Two functions are provided for each vector operation: a kernel for 
* performing the operation on the GPU, and a wrapper to be called by
* the host that calls the kernel internally. The kernel and its 
* wrapper have the same name, but the kernel name is preceded by an
* underscore (e.g., _addVV and addVV are the kernel and wrapper, 
* respectively). The function arguments appear in the same order in
* the wrapper and the kernel, but vectors are provided to the kernel
* as a bare pointer while they are provided to the wrapper as a 
* DeviceDArray object (passed by reference). The kernels require one
* input parameter that is not needed by the wrapper, namely the 
* vector length. 
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
* Vector assignment, a[i] = b[i], GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void eqV(DeviceDArray<cudaReal>& a, 
                  DeviceDArray<cudaReal> const & b);

/**
* Vector assignment, a[i] = b[i], GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void eqV(DeviceDArray<cudaComplex>& a, 
                  DeviceDArray<cudaComplex> const & b);

/**
* Vector assignment, a[i] = b, GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void eqS(DeviceDArray<cudaReal>& a, cudaReal const b);

/**
* Vector assignment, a[i] = b, GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void eqS(DeviceDArray<cudaComplex>& a, cudaComplex const b);


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
* Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _addVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n);

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
* Vector addition, a[i] = b[i] + c, GPU kernel (mixed type, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel (mixed type, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void addVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c);

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c);

/**
* Vector addition, a[i] = b[i] + c[i], GPU kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void addVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, cudaComplex const c);

/**
* Vector addition, a[i] = b[i] + c, GPU kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c);


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
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed type, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _subVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed type, c = real).
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
* Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed type, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed type, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void subVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c);

/**
* Vector subtraction, a[i] = b[i] - c[i], GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c);

/**
* Vector subtraction, a[i]=b[i]-c[i], kernel wrapper (mixed type, b=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c);

/**
* Vector subtraction, a[i]=b[i]-c[i], kernel wrapper (mixed type, c=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void subVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, cudaComplex const c);

/**
* Vector subtraction, a[i] = b[i] - c, GPU kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c);


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
* Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param n  size of arrays
*/
__global__ void _mulVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n);

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
* Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed type, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _mulVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n);

/**
* Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed type, c = real).
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
*/
__host__ void mulVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c);

/**
* Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void mulVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, cudaComplex const c);

/**
* Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c);


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
* Vector division, a[i] = b[i] / c[i], GPU kernel (mixed type, c = real).
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
* Vector division, a[i] = b[i] / c, GPU kernel (mixed type, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n);

/**
* Vector division, a[i] = b[i] / c[i], GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void divVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c);

/**
* Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed type, c=real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
*/
__host__ void divVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c);

/**
* Vector division, a[i] = b[i] / c, GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void divVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c);

/**
* Vector division, a[i] = b[i] / c, GPU kernel wrapper (mixed, c = real).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void divVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, cudaReal const c);


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
* Vector exponentiation, a[i] = exp(b[i]), GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void expV(DeviceDArray<cudaReal>& a, 
                   DeviceDArray<cudaReal> const & b);

/**
* Vector exponentiation, a[i] = exp(b[i]), GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void expV(DeviceDArray<cudaComplex>& a, 
                   DeviceDArray<cudaComplex> const & b);


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
* Vector addition in-place, a[i] += b[i], GPU kernel (mixed type).
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
* Vector addition in-place, a[i] += b, GPU kernel (mixed type).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _addEqS(cudaComplex* a, cudaReal const b, const int n);

/**
* Vector addition in-place, a[i] += b[i], GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void addEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector addition in-place, a[i] += b[i], GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void addEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaComplex> const & b);

/**
* Vector addition in-place, a[i] += b[i], GPU kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void addEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector addition in-place, a[i] += b, GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void addEqS(DeviceDArray<cudaReal>& a, cudaReal const b);

/**
* Vector addition in-place, a[i] += b, GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void addEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b);

/**
* Vector addition in-place, a[i] += b, GPU kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void addEqS(DeviceDArray<cudaComplex>& a, cudaReal const b);


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
* Vector subtraction in-place, a[i] -= b[i], GPU kernel (mixed type).
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
* Vector subtraction in-place, a[i] -= b, GPU kernel (mixed type).
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
*/
__host__ void subEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void subEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaComplex> const & b);

/**
* Vector subtraction in-place, a[i] -= b[i], kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void subEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void subEqS(DeviceDArray<cudaReal>& a, cudaReal const b);

/**
* Vector subtraction in-place, a[i] -= b, GPU kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void subEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b);

/**
* Vector subtraction in-place, a[i] -= b, GPU kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void subEqS(DeviceDArray<cudaComplex>& a, cudaReal const b);


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
* Vector multiplication in-place, a[i] *= b[i], GPU kernel (mixed type).
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
* Vector multiplication in-place, a[i] *= b, GPU kernel (mixed type).
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
*/
__host__ void mulEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void mulEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaComplex> const & b);

/**
* Vector multiplication in-place, a[i]*=b[i], kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void mulEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void mulEqS(DeviceDArray<cudaReal>& a, cudaReal const b);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void mulEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b);

/**
* Vector multiplication in-place, a[i] *= b, kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void mulEqS(DeviceDArray<cudaComplex>& a, cudaReal const b);


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
* Vector division in-place, a[i] /= b[i], GPU kernel (mixed type).
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
* Vector division in-place, a[i] /= b, GPU kernel (mixed type).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divEqS(cudaComplex* a, cudaReal const b, const int n);

/**
* Vector division in-place, a[i] /= b[i], GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void divEqV(DeviceDArray<cudaReal>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector division in-place, a[i] /= b[i], GPU kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void divEqV(DeviceDArray<cudaComplex>& a, 
                     DeviceDArray<cudaReal> const & b);

/**
* Vector division in-place, a[i] /= b, GPU kernel wrapper (cudaReal).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void divEqS(DeviceDArray<cudaReal>& a, cudaReal const b);

/**
* Vector division in-place, a[i] /= b, GPU kernel wrapper (mixed type).
*
* \param a  output array (LHS)
* \param b  input scalar (RHS)
*/
__host__ void divEqS(DeviceDArray<cudaComplex>& a, cudaReal const b);


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
__global__ void _addEqMulVS(cudaReal* a, cudaReal const * b, 
                            cudaReal const c, const int n);

/**
* Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
__host__ void addEqMulVS(DeviceDArray<cudaReal>& a, 
                         DeviceDArray<cudaReal> const & b, 
                         cudaReal const c);

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
*/
__host__ void sqNormV(DeviceDArray<cudaReal>& a, 
                      DeviceDArray<cudaComplex> const & b);

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
* \param n  size of arrays
*/
__global__ void _expMulVS(cudaReal* a, cudaReal const * b, 
                          cudaReal const c, const int n);

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
__host__ void expMulVS(DeviceDArray<cudaReal>& a, 
                       DeviceDArray<cudaReal> const & b, 
                       cudaReal const c);

/** @} */

} // namespace Vec
} // namespace Pscf
#endif
