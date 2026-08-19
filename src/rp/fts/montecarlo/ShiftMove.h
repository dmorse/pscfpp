#ifndef RP_SHIFT_MOVE_H
#define RP_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/fts/montecarlo/ShiftMove_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/fts/montecarlo/ShiftMove_u.h>
#endif

#endif
