#ifndef PSCF_CUDA_ERROR_CHECK_H
#define PSCF_CUDA_ERROR_CHECK_H

#include <util/global.h>
#include <string>
#include <cuda_runtime.h>

namespace Pscf {

/**
* Check if CUDA API function returns cudaSuccess (preprocessor macro).
*
* CUDA functions like cudaMalloc and cudaMemcpy return a cudaError_t
* object. If the operation was successful, this cudaError_t object 
* will be equal to cudaSuccess, otherwise the object will have some
* other value that indicates the type of error that occurred. If an
* error has occurred, this function prints the error string and throws
* an Exception.
*
* Kernel calls do not return a cudaError_t value, but a CUDA API 
* function can be called immediately after a kernel call to ensure 
* that the kernel was successfully launched. This may look like
* \code
*    cudaErrorCheck( cudaGetLastError() );
* \endcode
* However, CUDA kernels are executed asynchronously from the host 
* code, meaning that the kernel may not have even begun when 
* cudaGetLastError() is called. Therefore, this tool can only debug
* errors in a kernel launch, not in kernel execution.
*
* If one wished to perform error checks that test whether a kernel
* successfully executed, the following line of code could be used:
* \code
*    cudaErrorCheck( cudaDeviceSynchronize() );
* \endcode
* The cudaDeviceSynchronize() command pauses all operations on the 
* host until the GPU has completed executing all kernels, and then
* returns a cudaError_t object that indicates the error status of
* the device. These error checks are useful for debugging, but are 
* not recommended in most code because device synchronization can
* substantially reduce the overall speed of the code.
*
* \param ans  cudaError_t object returned from CUDA API function
*/
#define cudaErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__);}

/**
* Check if CUDA API function returns cudaSuccess (used by cudaErrorCheck macro).
*
* \param err  cudaError_t object returned from CUDA API function
* \param file  filename where cudaErrorCheck was called
* \param line  line number where cudaErrorCheck was called
*/
inline void gpuAssert(cudaError_t err, const char *file, int line)
{
   if (err != cudaSuccess) {
      std::string msg("CUDA Error: ");
      msg += cudaGetErrorString(err);
      throw Util::Exception(msg.c_str(), file, line);
   }
}

} // namespace Pscf
#endif
