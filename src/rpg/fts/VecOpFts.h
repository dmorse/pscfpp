#ifndef PSCF_VEC_OP_FTS_H
#define PSCF_VEC_OP_FTS_H

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/types.h>
#include <pscf/cuda/DeviceArray.h>

namespace Pscf {
namespace Rpg {

/** 
* Element-wise vector operations performed on the GPU for FTS classes.
*
* CUDA kernels that perform the operations are defined in VecOpFts.cu
* in an anonymous namespace, so they are not directly accessible. Kernel
* wrapper functions to be called by the host CPU, which call the kernel 
* internally, are public. 
*
* The output (the LHS of the vector operation) will always be the first
* parameter passed to the function. 
* 
* \ingroup Rpg_Fts_Module
* @{
*/
namespace VecOpFts {

   /**
   * Rescale array a from [0,1] to [-b, b], GPU kernel wrapper.
   *
   * \param a  array containing data between 0 and 1 (LHS)
   * \param b  scalar bound for rescaled array (RHS)
   */
   void mcftsScale(DeviceArray<Prdc::Cuda::cudaReal>& a, 
                   Prdc::Cuda::cudaReal const b);

   /**
   * Add array b to real part of a and array c to imaginary part of a
   *
   * \param a  output array of cudaComplex values
   * \param b  input array of cudaReal values
   * \param c  input array of cudaReal values
   */
   void fourierMove(DeviceArray<Prdc::Cuda::cudaComplex>& a, 
                    DeviceArray<Prdc::Cuda::cudaReal> const & b, 
                    DeviceArray<Prdc::Cuda::cudaReal> const & c);

   /**
   * Compute d field (functional derivative of H[w])
   */
   void computeDField(DeviceArray<Prdc::Cuda::cudaReal>& d, 
                      DeviceArray<Prdc::Cuda::cudaReal> const & Wc, 
                      DeviceArray<Prdc::Cuda::cudaReal> const & Cc, 
                      Prdc::Cuda::cudaReal const a, 
                      Prdc::Cuda::cudaReal const b, 
                      Prdc::Cuda::cudaReal const s);

   /**
   * Compute force bias
   */
   void computeForceBias(DeviceArray<Prdc::Cuda::cudaReal>& result, 
                         DeviceArray<Prdc::Cuda::cudaReal> const & di, 
                         DeviceArray<Prdc::Cuda::cudaReal> const & df, 
                         DeviceArray<Prdc::Cuda::cudaReal> const & dwc, 
                         Prdc::Cuda::cudaReal mobility);

}
/** @} */

}
}

#endif