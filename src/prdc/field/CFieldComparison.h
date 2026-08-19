#ifndef PRDC_C_FIELD_COMPARISON_H
#define PRDC_C_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// This is a unified header for a class template, which includes 
// the header for every partial specialization that is enabled
// for use by the build system.

#ifdef PSCF_CPP
#include <prdc/field/cpu/CFieldComparison.h>
#endif

#ifdef PSCF_CUDA
#include <prdc/field/cuda/CFieldComparison.h>
#endif

#endif
