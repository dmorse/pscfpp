#ifndef PSCF_VEC_OP_H
#define PSCF_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// This is a unified header for functions in the Pscf::VecOp namespace.
// It includes the header for every all functions that are compiled by
// the build system.

#ifdef PSCF_CPP
#include <pscf/backend/cpp/VecOpCx.h>
#endif

#ifdef PSCF_CUDA
#include <pscf/backend/cuda/VecOp.h>
#endif

#endif
