#ifndef PSCF_THREAD_ARRAY_H
#define PSCF_THREAD_ARRAY_H

#include <util/global.h>

namespace Pscf {

/**
* Global functions and variables to control GPU thread and block counts.
*/
namespace ThreadArray {

   /**
   * \defgroup Pscf_Cuda_ThreadArray_Module ThreadArray
   *
   * Management of GPU resources and setting of execution configurations.
   *
   * \ingroup Pscf_Cuda_Module
   * @{
   */
 
   /**
   * Initialize static variables in Pscf::ThreadArray namespace.
   */ 
   void init();

   /**
   * Set the number of threads per block to a default value.
   *
   * Query the hardware to determine a reasonable number.
   */
   void setThreadsPerBlock();

   /**
   * Set the number of threads per block to a specified value.
   * 
   * \param nThreadsPerBlock the number of threads per block (input)
   */
   void setThreadsPerBlock(int nThreadsPerBlock);

   /**
   * Given total number of threads, set 1D execution configuration.
   * 
   * Calculates the number of blocks and threads per block, and stores
   * the results in static variables along with the total number of 
   * threads requested. 
   * 
   * \param nThreadsLogical total number of required threads (input)
   */
   void setThreadsLogical(int nThreadsLogical);

   /**
   * Given total number of threads, set 1D execution configuration.
   *
   * Calculates the number of blocks and threads per block, and stores
   * the results in static variables along with the total number of 
   * threads requested. 
   * 
   * Sets the parameter nBlocks equal to the calculated number of blocks.
   * 
   * \param nThreadsLogical total number of required threads (input)
   * \param nBlocks  updated number of blocks (output)
   */
   void setThreadsLogical(int nThreadsLogical, int& nBlocks);

   /**
   * Given total number of threads, set 1D execution configuration.
   *
   * Calculates the number of blocks and threads per block, and stores
   * the results in static variables along with the total number of 
   * threads requested. 
   * 
   * Sets the parameter nBlocks equal to the calculated number of blocks,
   * and nThreads to the calculated number of threads per block.
   * 
   * \param nThreadsLogical total number of required threads (input)
   * \param nBlocks  updated number of blocks (output)
   * \param nThreads  updated number threads per block (output)
   */
   void setThreadsLogical(int nThreadsLogical, int& nBlocks, int& nThreads);

   /**
   * Check the execution configuration (threads and block counts).
   *
   * Check for validity and optimality, based on hardware warp size and 
   * streaming multiprocessor constraints. 
   */
   void checkExecutionConfig();

   // Accessors

   /**
   * Get the current number of blocks for execution.
   */
   int nBlocks();

   /**
   * Get the number of threads per block for execution.
   */
   int nThreads();

   /**
   * Return previously requested total number of threads.
   */
   int nThreadsLogical();

   /**
   * Get the warp size.
   */
   int warpSize();

   /**
   * Indicates whether there will be unused threads. 
   *
   * Returns true iff nThreads*nBlocks != nThreadsLogical.
   */
   bool hasUnusedThreads();

   /** @} */

} // namespace ThreadArray
} // namespace Pscf
#endif
