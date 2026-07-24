#ifndef RP_MIXTURE_H
#define RP_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/solvers/Mixture_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/solvers/Mixture_u.h>
#endif

#endif
