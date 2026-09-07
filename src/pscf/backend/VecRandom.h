#ifndef PSCF_VEC_RANDOM_H
#define PSCF_VEC_RANDOM_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// This is a unified header for a class template, which includes 
// the header for every explicit specialization that is enabled
// for use by the build system.

#ifdef PSCF_CPP
#include <pscf/backend/cpp/VecRandom.h>
#endif

#ifdef PSCF_CUDA
#include <pscf/backend/cuda/VecRandom.h>
#endif

#endif
