#ifndef RP_BINARY_STRUCTURE_FACTOR_H
#define RP_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef PSCF_CPP
#include <rp/fts/analyzer/BinaryStructureFactor_c.h>
#endif

#ifdef PSCF_CUDA
#include <rp/fts/analyzer/BinaryStructureFactor_u.h>
#endif

#endif
