#ifndef PSCF_THREADMESH_H
#define PSCF_THREADMESH_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <pscf/math/IntVec.h>
#include <vector_types.h> // from the CUDA library, defines dim3 type

namespace Pscf {

/**
* Management of multidimensional GPU thread execution configurations.
*
* For calculations involving data defined on a multidimensional mesh 
* (e.g. a field on a real-space grid), it is often useful to use a 
* multidimensional grid of CUDA threads defined on that same mesh
* (or on a slightly larger mesh if the mesh dimensions are not nicely
* divisible by powers of 2). This way, within the CUDA kernel, it is 
* easy to identify the grid point that corresponds to a given thread 
* using threadIdx and blockIdx. 
*/
namespace ThreadMesh {

   /**
   * \defgroup Pscf_Cuda_ThreadMesh_Module ThreadMesh
   *
   * Management of multidimensional GPU thread execution configurations.
   *
   * \ingroup Pscf_Cuda_Module
   * @{
   */

   /**
   * Given a multidimensional grid of threads, set execution configuration.
   * 
   * This function accepts an IntVec<D> array containing D mesh dimensions,
   * and constructs an execution configuration that is a D-dimensional grid 
   * of threads matching this mesh.
   * 
   * The function first determines an optimal number of threads per block,
   * based on maximizing thread occupancy. It then constructs a D-dimensional
   * thread block of that size (or less, if the requested grid of threads
   * is smaller than the optimal block size). Finally, it determines the 
   * dimensions of the D-dimensional grid of blocks that is needed to span
   * the entire input mesh.
   * 
   * The resulting grid and block dimensions are stored internally to the
   * ThreadMesh namespace, and can be accessed via the accessors gridDims() 
   * and blockDims, respectively. The dimensions are stored in objects of
   * type dim3 because that is the data type that is needed to launch a
   * CUDA kernel with a multidimensional mesh of threads. If D = 2, the
   * z member of these dim3 objects is set to 1, and if D = 1, the y 
   * member is also set to 1. 
   * 
   * IMPORTANT: A required input for this function is the boolean "invert".
   * If invert == true, the order of the dimensions stored in the dim3 
   * objects will be inverted from the dimensions that were provided in
   * meshDims. This means that
   *  - in 3D, the x, y, and z dimensions in meshDims will correspond to 
   *    the z, y, and x members of the dim3 objects, respectively.
   *  - in 2D, the x and y dimensions in meshDims will correspond to the
   *    y and x members of the dim3 objects, respectively.
   *  - in 1D, the result will be the same whether invert is true or false.
   * The reason that this option is provided is to more easily facilitate
   * coalesced data accessing in the CUDA kernel, explained further below.
   * 
   * If each thread in a warp (32 threads) must access a unique element of 
   * an array, that corresponds to 32 independent memory accesses. However,
   * if the elements being accessed by the warp are stored in consecutive
   * locations in memory, the compiler can "coalesce" this operation into 
   * one memory access, which can greatly speed up a CUDA kernel. 
   * 
   * A multidimensional thread block must be linearized into a 1D array 
   * before warps can be assigned. This linearization is performed in 
   * column-major order, in which x is the most rapidly changing index,
   * then y, then z. If, say, a thread block has 32 threads along it's x
   * dimension, then each warp will be at a single value of y and z, but
   * will have x values that span from 0 to 31.
   * 
   * Multidimensional arrays of data stored on the device must also be
   * linearized, but may be linearized in row-major order, in which z is
   * the most rapidly changing index, then y, then x. If the thread grid
   * for the CUDA kernel has the same dimensions as an array that was 
   * linearized in row-major order, the memory accesses within a single
   * warp cannot be coalesced. 
   * 
   * Instead, one can simply invert the dimensions within the CUDA kernel.
   * If the x dimension within the kernel corresponds to the z dimension
   * within the array, then an array that was linearized in row-major 
   * order can be accessed with coalescence. The CUDA kernel must be 
   * written with the knowledge that the dimensions are inverted in such
   * a way, but otherwise no changes are necessary.
   * 
   * Regardless of the value of "invert", the program will prioritize
   * having a block size of at least 32 in the x dimension if possible,
   * allowing coalescence to be maximized within a kernel.
   * 
   * \param meshDims  dimensions of the desired grid of threads (input)
   * \param invert  should the coordinate order be inverted, xyz -> zyx?
   * \param blockSize  desired block size (optional, must be power of 2)
   */
   template <int D>
   void setConfig(IntVec<D> const & meshDims, bool invert, 
                  int blockSize = -1);

   /**
   * Check the execution configuration (grid and block dimensions).
   *
   * Check for validity and optimality, based on hardware warp size and 
   * streaming multiprocessor constraints. 
   */
   void checkConfig();

   /**
   * Manually set the block size that should be used by default.
   * 
   * \param blockSize  the block size to be used
   */
   void setThreadsPerBlock(int blockSize);

   // Accessors

   /**
   * Get the multidimensional grid of blocks determined by setConfig.
   */
   dim3 gridDims();

   /**
   * Get the dimensions of a single block determined by setConfig.
   */
   dim3 blockDims();

   /**
   * Return last requested multidimensional grid of threads.
   */
   dim3 meshDims();

   /**
   * Get the warp size.
   */
   int warpSize();

   /**
   * Will there be unused threads?
   */
   bool hasUnusedThreads();

   /** @} */

} // namespace ThreadMesh
} // namespace Pscf
#endif