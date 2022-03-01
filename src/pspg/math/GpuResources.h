#ifndef PSPG_GPU_RESOURCES_H
#define PSPG_GPU_RESOURCES_H

#include "GpuTypes.h"            // typedefs used in Pspg
#include "ThreadGrid.h"          // Functions and variables for management of GPU execution configuration
#include "LinearAlgebra.h"       // Linear algebra kernels used in Pspg
#include "ParallelReductions.h"  // Kernels using parallel reduction algorithms
#include "KernelWrappers.h"      // Host functions in Pscf::Pspg namespace that wrap certain tricky kernels
#include "Helpers.h"             // Kernels very specific to certain calculations

#endif