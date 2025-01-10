#ifndef PSCF_VEC_OP_FTS_H
#define PSCF_VEC_OP_FTS_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/GpuResources.h>

namespace Pscf {
namespace Rpg {
namespace VecOpFts {

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

   /**
   * Rescale array a from [0,1] to [-b, b], GPU kernel wrapper.
   *
   * \param a  array containing data between 0 and 1 (LHS)
   * \param b  scalar bound for rescaled array (RHS)
   */
   void mcftsScale(DeviceArray<cudaReal>& a, cudaReal const b);

   /**
   * Add array b to real part of a and array c to imaginary part of a
   *
   * \param a  output array of cudaComplex values
   * \param b  input array of cudaReal values
   * \param c  input array of cudaReal values
   */
   void fourierMove(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c);

   /**
   * Compute d field (functional derivative of H[w])
   */
   void computeDField(DeviceArray<cudaReal>& d, 
                      DeviceArray<cudaReal> const & Wc, 
                      DeviceArray<cudaReal> const & Cc, 
                      cudaReal const a, cudaReal const b, cudaReal const s);

   /**
   * Compute force bias
   */
   void computeForceBias(DeviceArray<cudaReal>& result, 
                         DeviceArray<cudaReal> const & di, 
                         DeviceArray<cudaReal> const & df, 
                         DeviceArray<cudaReal> const & dwc, cudaReal mobility);

   /** @} */

}
}
}

#endif