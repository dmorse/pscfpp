#ifndef PSCF_CUDA_RESOURCES_H
#define PSCF_CUDA_RESOURCES_H

#include <pscf/backend/cuda/VecOp.h>        // Vector operations 
#include <pscf/backend/cuda/Reduce.h>       // Parallel reduction algorithms
#include <pscf/backend/cuda/cudaTypes.h>    // typedefs for real and complex 
#include <pscf/backend/cuda/DeviceArray.h>  // Array container stored on device
#include <pscf/backend/cuda/HostDArray.h>   // Array container stored on host
#include <pscf/backend/cuda/ThreadArray.h>  // Manager of GPU execution configuration
#include <pscf/backend/cuda/ThreadMesh.h>   // Manager of GPU execution configuration

#endif
