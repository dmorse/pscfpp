#ifndef PSCF_VEC_OP_H
#define PSCF_VEC_OP_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

// This is a unified header for functions in the Pscf::VecOp namespace.
// It includes the header for every all functions that are compiled by
// the build system.

#include <pscf/backend/cpp/VecOp.h>

#ifdef PSCF_CUDA
#include <pscf/backend/cuda/VecOp.h>
#endif

#endif
