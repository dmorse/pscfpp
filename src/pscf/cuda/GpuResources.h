#ifndef PSCF_GPU_RESOURCES_H
#define PSCF_GPU_RESOURCES_H

#include "GpuTypes.h"           // typedefs used in cuda code
#include "ThreadGrid.h"         // Management of GPU execution configuration
#include "LinearAlgebra.h"      // Linear algebra kernels 
#include "VecOp.h"              // Vector operation kernels and wrappers
#include "ParallelReductions.h" // Kernels using parallel reduction algorithms
#include "KernelWrappers.h"     // Host functions for reductions
#include "DeviceArray.h"        // Array container stored on device
#include "HostDArray.h"         // Array container stored on host

#endif
