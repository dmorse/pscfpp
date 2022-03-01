#ifndef PSPG_THREADGRID_H
#define PSPG_THREADGRID_H

#include "GpuTypes.h"
#include <util/global.h>

namespace Pscf {
namespace Pspg {
namespace ThreadGrid {

  /**
  * \defgroup Pspg_ThreadGrid_Module ThreadGrid
  *
  * Management of GPU resources and setting of execution configurations.
  *
  * \ingroup Pscf_Pspg_Module
  * @{
  */
  
  void init();

  /**
  * Set the number of threads per block by querying the 
  * hardware to determine a reasonable number.
  */
  void setThreadsPerBlock();

  /**
  * Set the number of threads per block via input from the
  * program or user.
  * 
  * \param nThreadsPerBlock the number of threads per block 
  */
  void setThreadsPerBlock(int const nThreadsPerBlock);

  /**
  * Set the total number of threads required for execution. Recalculate the
  * number of blocks and threads per block if necessary.
  * 
  * \param nThreadsLogical the total number of threads required for execution
  */
  void setThreadsLogical(int const nThreadsLogical);

  /**
  * Set the total number of threads required for execution. Recalculate the
  * number of blocks and threads per block if necessary. Also updates the 
  * inputted nBlocks parameter with the corresponding value.
  * 
  * \param nThreadsLogical the total number of threads required for execution
  * \param nBlocks passed by reference, updated with the number of blocks for execution
  */
  void setThreadsLogical(int const nThreadsLogical, int & nBlocks);

  /**
  * Set the total number of threads required for execution. Recalculate the
  * number of blocks and threads per block if necessary. Also updates the 
  * inputted nBlocks and nThreads parameters with the corresponding values.
  * 
  * \param nThreadsLogical the total number of threads required for execution
  * \param nBlocks passed by reference, updated with the number of blocks for execution
  * \param nThreads passed by referenced, updated with the threads per block for execution
  */
  void setThreadsLogical(int const nThreadsLogical, int & nBlocks, int & nThreads);

  /**
  * Check the execution configuration (number of threads and threads per block) for validity 
  * and optimality, based on hardware warp size and streaming multiprocessor constraints. 
  */
  void checkExecutionConfig();

  // Accessors

  /// Get the number of blocks for execution.
  int nBlocks();

  /// Get the number of threads per block for execution.
  int nThreads();

  /// Return previously requested total number of threads.
  int nThreadsLogical();

  /// Indicates whether there will be unused threads. 
  /// This is the case if nThreads*nBlocks != nThreadsLogical.
  bool hasUnusedThreads();

  /** @} */

}
}
}
#endif