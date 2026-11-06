#ifndef PSCF_HOST_SEND_ARRAY_H
#define PSCF_HOST_SEND_ARRAY_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

// This is a unified header for a class template, which includes 
// the header for every explicit specialization that is enabled
// for use by the build system.

#include <pscf/backend/cpp/HostSendArray.h>

#ifdef PSCF_CUDA
#include <pscf/backend/cuda/HostSendArray.h>
#endif

#endif
