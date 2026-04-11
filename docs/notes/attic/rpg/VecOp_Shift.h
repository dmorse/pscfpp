#ifndef PSCF_VEC_OP_FTS_H
#define PSCF_VEC_OP_FTS_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/math/IntVec.h>
#include <vector_types.h>

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
   * Shift w Field
   */
   template <int D>
   void shiftWField(DeviceArray<cudaReal>& wshift,
                    DeviceArray<cudaReal>const & w0,
                    IntVec<D> const & meshDims,
                    IntVec<D> shift);

}
/** @} */

}
}

#endif
