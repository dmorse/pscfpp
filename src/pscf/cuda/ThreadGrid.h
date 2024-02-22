#ifndef PSCF_THREADGRID_H
#define PSCF_THREADGRID_H

#include "GpuTypes.h"
#include <util/global.h>

namespace Pscf {

/**
* Global functions and variables to control GPU thread and block counts.
*/
namespace ThreadGrid {

  /**
  * \defgroup Pscf_Cuda_ThreadGrid_Module ThreadGrid
  *
  * Management of GPU resources and setting of execution configurations.
  *
  * \ingroup Pscf_Cuda_Module
  * @{
  */
 
  /**
  * Initialize static variables in Pscf::ThreadGrid namespace.
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
  * Set the total number of threads required for execution. 
  * 
  * Calculate the number of blocks, and calculate threads per block if 
  * necessary. Updates static variables.
  * 
  * \param nThreadsLogical total number of required threads (input)
  */
  void setThreadsLogical(int nThreadsLogical);

  /**
  * Set the total number of threads required for execution. 
  *
  * Recalculate the number of blocks, and calculate threads per block if 
  * necessary. Also updates the nBlocks output parameter.
  * 
  * \param nThreadsLogical total number of required threads (input)
  * \param nBlocks  updated number of blocks (output)
  */
  void setThreadsLogical(int nThreadsLogical, int& nBlocks);

  /**
  * Set the total number of threads required for execution. 
  *
  * Computes and sets the number of blocks, and sets threads per block 
  * if necessary.  Updates values of nBlocks and nThreads parameters in
  * output parameters that are passed by value.
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
  * Indicates whether there will be unused threads. 
  *
  * Returns true iff nThreads*nBlocks != nThreadsLogical.
  */
  bool hasUnusedThreads();

  /** @} */

} // namespace ThreadGrid
} // namespace Pscf
#endif
