#ifndef RP_BLOCK_H
#define RP_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/solvers/Block_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/solvers/Block_u.h>
#endif

#endif
