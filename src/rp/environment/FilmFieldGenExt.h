#ifndef RP_FILM_FIELD_GEN_EXT_H
#define RP_FILM_FIELD_GEN_EXT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/environment/FilmFieldGenExt_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/environment/FilmFieldGenExt_u.h>
#endif

#endif
