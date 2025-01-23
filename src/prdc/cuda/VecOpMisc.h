#ifndef PRDC_VEC_OP_MISC_H
#define PRDC_VEC_OP_MISC_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "types.h"
#include <pscf/cuda/DeviceArray.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {
namespace VecOp {

/*
* Miscellaneous element-wise vector operations performed on the GPU.
*
* Note: this file is included at the end of VecOp.h, so any file that 
* includes VecOp.h will also include this file.
*
* The functions defined in this file are wrappers for CUDA kernels that 
* perform the actual vector operations. The kernels themselves are only 
* intended to be called through their wrappers, so they are defined in 
* an anonymous namespace in VecOpMisc.cu.
*
* The functions defined in this file combine 2 or more element-wise 
* vector operations into a single kernel launch, which will perform the 
* operation faster than consecutively calling multiple of the functions 
* in VecOp.h. These functions are not intended to be comprehensive. 
* Rather, they are written and included as needed during the development 
* of other code.
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
*/

// Functions that combine multiple VecOp operations
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Vector addition w/ coefficient, a[i] = (b[i]*c) + (d[i]*e), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array 1 (RHS)
* \param c  input scalar (RHS)
* \param d  input array 2 (RHS)
* \param e  input scalar 2 (RHS)
*/
void addVcVc(DeviceArray<cudaReal>& a, 
             DeviceArray<cudaReal> const & b, cudaReal const c, 
             DeviceArray<cudaReal> const & d, cudaReal const e);

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
void addVcVcVc(DeviceArray<cudaReal>& a, 
               DeviceArray<cudaReal> const & b, cudaReal const c, 
               DeviceArray<cudaReal> const & d, cudaReal const e, 
               DeviceArray<cudaReal> const & f, cudaReal const g);

/**
* Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
void addEqVc(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
             cudaReal const c);

/**
* Vector subtraction, a[i] = b[i] - c[i] - d, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input array (RHS)
* \param d  input scalar (RHS)
*/
void subVVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
            DeviceArray<cudaReal> const & c, cudaReal const d);

/**
* Vector division in-place w/ coeff., a[i] /= (b[i] * c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar (RHS)
*/
void divEqVc(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
             cudaReal const c);

/**
* Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
* \param c  input scalar
*/
void expVc(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           cudaReal const c);


// Pair functions
// ~~~~~~~~~~~~~~

/**
* Vector assignment in pairs, ax[i] = b[i], kernel wrapper.
*
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param s  shared input array to be assigned to both a1 and a2
*/
void eqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
             DeviceArray<cudaReal> const & s);

/**
* Vector multiplication in pairs, ax[i] = bx[i] * s[i], kernel wrapper.
* 
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param b1  input array 1 (RHS)
* \param b2  input array 2 (RHS)
* \param s  shared input array to be multiplied by both b1 and b2
*/
void mulVVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
               DeviceArray<cudaReal> const & b1, 
               DeviceArray<cudaReal> const & b2, 
               DeviceArray<cudaReal> const & s);

/**
* In-place vector multiplication in pairs, ax[i] *= s[i], kernel wrapper
* 
* \param a1  output array 1 (LHS)
* \param a2  output array 2 (LHS)
* \param s  shared input array to be multiplied by both a1 and a2
*/
void mulEqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
                DeviceArray<cudaReal> const & s);


// Functions of "many" vectors
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Add an undefined number of vectors pointwise, kernel wrapper.
*
* The input array 'vecs' contains the arrays that will be added together. 
* The size of vecs determines the number of vectors that will ultimately
* be added together by the GPU kernel.
*
* \param a  output array (LHS)
* \param vecs  array of DeviceArrays to be added
*/
void addVMany(DeviceArray<cudaReal>& a, 
              DArray<DeviceArray<cudaReal> > const & vecs);

/**
* Add an undefined number of vectors pointwise, kernel wrapper.
*
* The input array 'vecs' contains const pointers to each array that will be 
* added together. The size of vecs determines the number of vectors that
* will ultimately be added together by the GPU kernel.
*
* This version of addVMany is provided for cases in which one needs to add
* many arrays that are not already stored together in a DArray. The caller 
* must simply assemble an array of pointers to all of the arrays that should
* be added, and then pass it to this method.
*
* \param a  output array (LHS)
* \param vecs  array of pointers to DeviceArrays to be added
*/
void addVMany(DeviceArray<cudaReal>& a, 
              DArray<DeviceArray<cudaReal> const *> const & vecs);

/**
* Multiply an undefined number of vectors pointwise, kernel wrapper.
*
* The input array 'vecs' contains the arrays that will be multiplied. 
* The size of vecs determines the number of vectors that will ultimately
* be multiplied together by the GPU kernel.
*
* \param a  output array (LHS)
* \param vecs  array of DeviceArrays to be multiplied
*/
void mulVMany(DeviceArray<cudaReal>& a, 
              DArray<DeviceArray<cudaReal> > const & vecs);

/**
* Multiply an undefined number of vectors pointwise, kernel wrapper.
*
* The input array 'vecs' contains const pointers to each array that will be 
* multiplied together. The size of vecs determines the number of vectors
* that will ultimately be multiplied together by the GPU kernel.
* 
* This version of mulVMany is provided for cases in which one needs to 
* multiply many arrays that are not already stored together in a DArray. 
* The caller must simply assemble an array of pointers to all of the 
* arrays that should be multiplied, and then pass it to this method.
*
* \param a  output array (LHS)
* \param vecs  array of pointers to DeviceArrays to be multiplied
*/
void mulVMany(DeviceArray<cudaReal>& a, 
              DArray<DeviceArray<cudaReal> const *> const & vecs);


// Other useful functions
// ~~~~~~~~~~~~~~~~~~~~~~

/**
* Squared norm of complex number, a[i] = norm(b[i])^2, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
void sqNormV(DeviceArray<cudaReal>& a, DeviceArray<cudaComplex> const & b);

/**
* Norm of complex number to the 4th power, a[i] = norm(b[i])^4, kernel wrapper.
*
* \param a  output array (LHS)
* \param b  input array (RHS)
*/
void sqSqNormV(DeviceArray<cudaReal>& a, DeviceArray<cudaComplex> const & b);

/** @} */

} // namespace VecOp
} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif
