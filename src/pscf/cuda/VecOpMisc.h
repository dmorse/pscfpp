#ifndef PSCF_CUDA_VEC_OP_MISC_H
#define PSCF_CUDA_VEC_OP_MISC_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace VecOp {

   /*
   * Miscellaneous element-wise vector operations performed on the GPU.
   *
   * This collection of functions is not intended to be comprehensive.  
   * Rather, they are written and included as needed during the development 
   * of other code.
   *
   * Note: this file is included at the end of VecOp.h, so any file that
   * includes VecOp.h will also include this file.
   *
   * Combined operations:
   *
   * The first set of functions defined in this file combine 2 or more 
   * element-wise vector operations. Several of these functions perform
   * linear combinations of vectors multiplied (or scaled) by coefficients, 
   * thus combining vector-scalar multiplication with vector addition. On 
   * a GPU, these operations are performed by launching a single kernel,
   * which will generally be faster than consecutively calling multiple 
   * simpler operations.
   *
   * The names of these functions follow conventions similar to those used
   * in VecOp.h, using eq, add, sub, mul, div, exp and combinations
   * thereof to indicate the operation(s) being performed. V denotes a
   * a vector input, S denotes a scalar, and Vc denotes a vector that 
   * is multiplied by a scalar coefficient and then used in another
   * operation. For example, addEqVc(a, b, c) performs a[i] += b[i] * c
   * for all i.
   *
   * Pair operations:
   * 
   * A second of functions defined in this file contain the word Pair,
   * indicating that these functions perform the same operation for a
   * pair of real output arrays, using a shared input array. For example, 
   * eqVPair performs a1[i] = s[i] and a2[i] = s[i] for all i. On a GPU,
   * performing these operations in pairs is faster because the shared 
   * input array only needs to be loaded from global memory once.
   *
   * "Many" operations:
   * 
   * A third set of functions defined in this file contain the word
   * "Many", indicating that an undefined number of input vectors (>2) 
   * are involved in an operation. For example, addVMany adds >2 vectors
   * together by passing an array of vectors, rather than a discrete set 
   * of vectors.
   *
   * The functions declared in this file are wrappers for CUDA kernels
   * that perform the actual vector operations. The underlying kernels 
   * are only intended to be called through their wrappers, and so are
   * defined in an anonymous namespace in the file VecOpMisc.cu that
   * is inaccessible outside that file. 
   */

   // Combined operations

   /**
   * Add two scaled vectors, a[i] = b1[i]*c1 + b2[i]*c2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b1  real input array 1 (RHS)
   * \param c1  real coefficent of b1 (RHS)
   * \param b2  real input array 2 (RHS)
   * \param c2  real coefficent of b2 (RHS)
   */
   void addVcVc(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b1, cudaReal const c1,
                DeviceArray<cudaReal> const & b2, cudaReal const c2);

   /**
   * Add a scaled vector and a scalar, a[i] = b[i]*c + s (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real input array (RHS)
   * \param c  real coefficient of b (RHS)
   * \param s  real scalar summand (RHS)
   */
   void addVcS(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b, cudaReal const c,
               cudaReal const s);

   /**
   * Add a scaled vector in-place, a[i] += b[i] * c (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real output array (LHS)
   * \param b  real input array (RHS)
   * \param c  real scalar coefficient of b (RHS)
   */
   void addEqVc(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b,
                cudaReal const c);

   /**
   * Add 3 scaled vectors, a[i] = b1[i]*c1 + b2[i]*c2 + b3[i]*c3 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real output array (LHS)
   * \param b1  real input array 1 (RHS)
   * \param c1  real coefficient of b1 (RHS)
   * \param b2  real input array 2 (RHS)
   * \param c2  real coefficient of b2 (RHS)
   * \param b3  real input array 3 (RHS)
   * \param c3  real coefficient of b3 (RHS)
   */
   void addVcVcVc(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b1, cudaReal const c1,
                  DeviceArray<cudaReal> const & b2, cudaReal const c2,
                  DeviceArray<cudaReal> const & b3, cudaReal const c3);

   /**
   * Add 2 scaled vectors + scalar, a[i] = b1[i]*c1 + b2[i]*c2 + s (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real output array (LHS)
   * \param b1  real input array 1 (RHS)
   * \param c1  real coefficient of b1 (RHS)
   * \param b2  real input array 2 (RHS)
   * \param c2  real coefficient of b2 (RHS)
   * \param s  real scalar summand (RHS)
   */
   void addVcVcS(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaReal> const & b1, cudaReal const c1,
                 DeviceArray<cudaReal> const & b2, cudaReal const c2,
                 cudaReal const s);

   /**
   * Vector division in-place w/ coeff., a[i] /= (b[i] * c).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  complex array (LHS)
   * \param b  real array (RHS)
   * \param c  input scalar (RHS)
   */
   void divEqVc(DeviceArray<cudaComplex>& a,
                DeviceArray<cudaReal> const & b,
                cudaReal const c);

   /**
   * Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  real array (RHS)
   * \param c  input scalar
   */
   void expVc(DeviceArray<cudaReal>& a, 
              DeviceArray<cudaReal> const & b,
              cudaReal const c);


   // Pair operations (two output arrays and a shared input)

   /**
   * Vector assignment in pairs, ax[i] = b[i], x = 1, 2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a1  real array 1 (LHS)
   * \param a2  real array 2 (LHS)
   * \param s  shared real array to be assigned to both a1 and a2
   */
   void eqVPair(DeviceArray<cudaReal>& a1, 
                DeviceArray<cudaReal>& a2,
                DeviceArray<cudaReal> const & s);

   /**
   * Vector multiplication in pairs, ax[i] = bx[i] * s[i], x=1,2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a1  real array 1 (LHS)
   * \param a2  real array 2 (LHS)
   * \param b1  real array 1 (RHS)
   * \param b2  real array 2 (RHS)
   * \param s  shared real array to be multiplied by both b1 and b2
   */
   void mulVVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2,
                  DeviceArray<cudaReal> const & b1,
                  DeviceArray<cudaReal> const & b2,
                  DeviceArray<cudaReal> const & s);

   /**
   * In-place vector multiplication in pairs, ax[i] *= s[i], x=1,2 (real).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a1  real array 1 (LHS)
   * \param a2  real array 2 (LHS)
   * \param s  shared real array to multiply both a1 and a2 (RHS)
   */
   void mulEqVPair(DeviceArray<cudaReal>& a1, 
                   DeviceArray<cudaReal>& a2,
                   DeviceArray<cudaReal> const & s);

   // Functions of "many" vectors

   /**
   * Add an arbitrary number of vectors pointwise (real).
   *
   * The input array 'vecs' contains the arrays that will be added.
   * The size of array vecs determines the number of vectors that
   * will be added together.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of DeviceArrays to be added (RHS)
   */
   void addVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> > const & vecs);

   /**
   * Add an arbitrary number of vectors pointwise (real).
   *
   * The input array 'vecs' contains const pointers to each array that
   * will be added together. The size of vecs determines the number of
   * vectors that will be added together.
   *
   * This version of addVMany is provided for cases in which one needs to
   * add many arrays that are not already stored together in a DArray. The
   * caller must simply assemble an array of pointers to all of the arrays
   * that should be added, and then pass it to this method.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of pointers to DeviceArrays to be added
   */
   void addVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> const *> const & vecs);

   /**
   * Multiply an undefined number of vectors pointwise (real).
   *
   * The input array 'vecs' contains the arrays that will be multiplied.
   * The size of vecs determines the number of vectors that will be
   * multiplied together.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of DeviceArrays to be multiplied (RHS)
   */
   void mulVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> > const & vecs);

   /**
   * Multiply an undefined number of vectors pointwise.
   *
   * The input array 'vecs' contains const pointers to each array that
   * will be multiplied together. The size of vecs determines the number
   * of vectors that will ultimately be multiplied together pointwise.
   *
   * This version of mulVMany is provided for cases in which one needs to
   * multiply many arrays that are not already stored together in a DArray.
   * The caller must simply assemble an array of pointers to all of the
   * arrays that should be multiplied, and then pass it to this method.
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param vecs  array of real arrays to be multiplied (RHS)
   */
   void mulVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> const *> const & vecs);

   // Other useful functions

   /**
   * Fourth power of magnitude, a[i] = |b[i]|^4 (complex).
   *
   * \ingroup Pscf_Cuda_VecOp_Module
   *
   * \param a  real array (LHS)
   * \param b  complex array (RHS)
   */
   void sqSqAbsV(DeviceArray<cudaReal>& a,
                 DeviceArray<cudaComplex> const & b);

} // namespace VecOp
} // namespace Pscf
#endif
