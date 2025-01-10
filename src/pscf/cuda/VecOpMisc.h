#ifndef PSCF_VEC_OP_MISC_H
#define PSCF_VEC_OP_MISC_H

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
* Miscellaneous element-wise vector operations performed on the GPU.
*
* Note: this file is included at the end of VecOp.h, so any file that 
* includes VecOp.h will also include this file.
*
* CUDA kernels that perform the vector operations are defined here, as 
* well as wrapper functions to be called by the host CPU, which call 
* the kernel internally. The kernels should always be accessed through 
* the wrapper functions when possible. 
*
* These functions combine 2 or more element-wise vector operations into 
* a single kernel launch, which will perform the operation much faster
* than consecutively calling multiple of the functions in VecOp.h. These
* functions are not intended to be comprehensive. Rather, they are 
* written and included as needed during the development of other code.
*
* The names of these functions follow the same conventions as those in
* VecOp, using add, sub, mul, div, exp, eq, and combinations thereof to
* indicate the operation(s) being performed. V denotes a vector, S 
* denotes a scalar, and Vc denotes a vector that is multiplied by a 
* scalar coefficient and then used in another operation. For example, 
* addEqVc(a, b, c) performs a[i] += b[i] * c for all i. 
*
* Another set of functions defined in this file contain the word Pair, 
* indicating that these functions perform the same operation for a pair
* of output arrays. For example, eqVPair performs a1[i] = s[i] and 
* a2[i] = s[i] for all i. Performing these operations in pairs is 
* faster because the array s only needs to be loaded from global memory
* once. 
* 
* A third set of functions defined in this file contain the word "Many",
* indicating that an undefined number of vectors (>2) are involved in an 
* operation. For example, addVMany adds >2 vectors together by passing 
* an array of vectors, rather than a discrete set of vectors.
* 
* \ingroup Pscf_Cuda_Module 
* @{
*/


// Functions that combine multiple VecOp operations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector addition w/ coefficient, a[i] = (b[i]*c) + (d[i]*e), GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array 1 (RHS)
* \param c  input scalar 1 (RHS)
* \param d  input array 2 (RHS)
* \param e  input scalar 2 (RHS)
* \param n  size of arrays
*/
__global__ void _addVcVc(cudaReal* a, cudaReal const * b, cudaReal const c, 
                         cudaReal const * d, cudaReal const e, const int n);

/**
* Vector addition w/ coefficient, a[i] = (b[i]*c) + (d[i]*e), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array 1 (RHS)
* \param c  input scalar (RHS)
* \param d  input array 2 (RHS)
* \param e  input scalar 2 (RHS)
*/
__host__ void addVcVc(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaReal> const & b, cudaReal const c, 
                      DeviceArray<cudaReal> const & d, cudaReal const e);

/**
* 3-vector addition w/ coeff, a[i] = (b[i]*c) + (d[i]*e) + (f[i]*g), GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array 1 (RHS)
* \param c  input scalar 1 (RHS)
* \param d  input array 2 (RHS)
* \param e  input scalar 2 (RHS)
* \param f  input array 3 (RHS)
* \param g  input scalar 3 (RHS)
* \param n  size of arrays
*/
__global__ void _addVcVcVc(cudaReal* a, cudaReal const * b, cudaReal const c, 
                           cudaReal const * d, cudaReal const e, 
                           cudaReal const * f, cudaReal const g, const int n);

/**
* 3-vector add. w/ coeff, a[i] = (b[i]*c) + (d[i]*e) + (f[i]*g), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array 1 (RHS)
* \param c  input scalar 1 (RHS)
* \param d  input array 2 (RHS)
* \param e  input scalar 2 (RHS)
* \param f  input array 3 (RHS)
* \param g  input scalar 3 (RHS)
*/
__host__ void addVcVcVc(DeviceArray<cudaReal>& a, 
                        DeviceArray<cudaReal> const & b, cudaReal const c, 
                        DeviceArray<cudaReal> const & d, cudaReal const e, 
                        DeviceArray<cudaReal> const & f, cudaReal const g);

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
*/
__host__ void addEqVc(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaReal> const & b, 
                      cudaReal const c);

/**
* Vector subtraction, a[i] = b[i] - c[i] - d, GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param d  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _subVVS(cudaReal* a, cudaReal const * b, 
                        cudaReal const * c, cudaReal const d, const int n);

/**
* Vector subtraction, a[i] = b[i] - c[i] - d, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param d  input scalar (RHS)
*/
__host__ void subVVS(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b, 
                     DeviceArray<cudaReal> const & c, cudaReal const d);

/**
* Vector division in-place w/ coeff., a[i] /= (b[i] * c), GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
* \param n  size of arrays
*/
__global__ void _divEqVc(cudaComplex* a, cudaReal const * b, 
                         cudaReal const c, const int n);

/**
* Vector division in-place w/ coeff., a[i] /= (b[i] * c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
__host__ void divEqVc(DeviceArray<cudaComplex>& a, 
                      DeviceArray<cudaReal> const & b, cudaReal const c);

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
*/
__host__ void expVc(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    cudaReal const c);


// Pair functions
// ~~~~~~~~~~~~~~

/**
* Vector assignment in pairs, a1[i] = b[i] and a2[i] = b[i], CUDA kernel.
*
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param s  shared input array to be assigned to both a1 and a2
* \param n  size of arrays
*/
__global__ void _eqVPair(cudaReal* a1, cudaReal* a2, 
                         cudaReal const * s, const int n);

/**
* Vector assignment in pairs, ax[i] = b[i], kernel wrapper.
*
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param s  shared input array to be assigned to both a1 and a2
*/
__host__ void eqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
                      DeviceArray<cudaReal> const & s);

/**
* Vec. mul. in pairs, ax[i] = bx[i] * s[i], CUDA kernel.
* 
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param b1  input array 1 (RHS)
* \param b2  input array 2 (RHS)
* \param s  shared input array to be multiplied by both a1 and a2
* \param n  size of arrays
*/
__global__ void _mulVVPair(cudaReal* a1, cudaReal * a2, 
                           cudaReal const * b1, cudaReal const * b2, 
                           cudaReal const * s, const int n);

/**
* Vec. mul. in pairs, ax[i] = bx[i] * s[i], kernel wrapper.
* 
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param b1  input array 1 (RHS)
* \param b2  input array 2 (RHS)
* \param s  shared input array to be multiplied by both b1 and b2
*/
__host__ void mulVVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
                        DeviceArray<cudaReal> const & b1, 
                        DeviceArray<cudaReal> const & b2, 
                        DeviceArray<cudaReal> const & s);

/**
* In-place vec. mul. in pairs, ax[i] *= s[i], CUDA kernel
* 
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param s  shared input array to be multiplied by both a1 and a2
* \param n  size of arrays
*/
__global__ void _mulEqVPair(cudaReal* a1, cudaReal* a2, 
                            cudaReal const * s, const int n);

/**
* In-place vec. mul. in pairs, ax[i] *= s[i], kernel wrapper
* 
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param s  shared input array to be multiplied by both a1 and a2
*/
__host__ void mulEqVPair(DeviceArray<cudaReal>& a1, 
                         DeviceArray<cudaReal>& a2, 
                         DeviceArray<cudaReal> const & s);


// Functions of "many" vectors
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


// Other useful functions
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Squared norm of complex number, a[i] = norm(b[i])^2, GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _sqNormV(cudaReal* a, cudaComplex const * b, const int n);

/**
* Squared norm of complex number, a[i] = norm(b[i])^2, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void sqNormV(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaComplex> const & b);

/**
* Norm of complex number to the 4th power, a[i] = norm(b[i])^4, GPU kernel.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param n  size of arrays
*/
__global__ void _sqSqNormV(cudaReal* a, cudaComplex const * b, const int n);

/**
* Norm of complex number to the 4th power, a[i] = norm(b[i])^4, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
__host__ void sqSqNormV(DeviceArray<cudaReal>& a, 
                        DeviceArray<cudaComplex> const & b);

/** @} */

} // namespace VecOp
} // namespace Pscf
#endif
