#ifndef PSCF_CUDA_ERROR_CHECK_H
#define PSCF_CUDA_ERROR_CHECK_H

#include <util/global.h>
#include <string>
#include <cuda_runtime.h>

namespace Pscf {

#define cudaErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__);}

inline void gpuAssert(cudaError_t code, const char *file, int line)
{
   if (code != cudaSuccess) {
      std::string msg("CUDA Error: ");
      msg += cudaGetErrorString(code);
      throw Util::Exception(msg.c_str(), file, line);
   }
}

} // namespace Pscf
#endif
