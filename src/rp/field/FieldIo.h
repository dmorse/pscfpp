#ifndef RP_FIELD_IO_H
#define RP_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/field/FieldIo_cp.h>
#endif

#ifdef PSCF_CUDA
#include <rp/field/FieldIo_cu.h>
#endif

#endif
