#ifndef RP_MAX_ORDER_PARAMETER_H
#define RP_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/fts/analyzer/MaxOrderParameter_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/fts/analyzer/MaxOrderParameter_u.h>
#endif

#endif
