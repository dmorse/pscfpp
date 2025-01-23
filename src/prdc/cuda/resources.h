#ifndef PSCF_GPU_RESOURCES_H
#define PSCF_GPU_RESOURCES_H

#include "types.h"         // typedefs used in cuda code
#include "VecOp.h"         // Vector operation kernels and wrappers
#include "Reduce.h"        // Functions using parallel reduction algorithms
#include <pscf/cuda/DeviceArray.h> // Array container stored on device
#include <pscf/cuda/HostDArray.h>  // Array container stored on host
#include <pscf/cuda/ThreadGrid.h>  // Management of GPU execution configuration

#endif
