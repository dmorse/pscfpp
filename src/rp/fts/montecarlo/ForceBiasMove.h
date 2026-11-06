#ifndef RP_FORCE_BIAS_MOVE_H
#define RP_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/fts/montecarlo/ForceBiasMove_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/fts/montecarlo/ForceBiasMove_u.h>
#endif

#endif
