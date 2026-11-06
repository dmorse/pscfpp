#ifndef PSCF_WAVE_LIST_H
#define PSCF_WAVE_LIST_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

// This is a unified header for a class template, which includes 
// the header for every partial specialization that is enabled
// for use by the build system.

#ifdef PSCF_CPP
#include <pscf/field/cpp/WaveList.h>
#endif

#ifdef PSCF_CUDA
#include <pscf/field/cuda/WaveList.h>
#endif

#endif
