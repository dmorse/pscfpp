#ifndef RPG_AM_COMP_BASE_H
#define RPG_AM_COMP_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* Explicit instantiation declarations for specializations of AmIteratorTmpl
* that are used as base classes for specializations of Rpg::AmCompressor
* and Rpg::LrAmCompressor.
*/

#include <pscf/iterator/AmIteratorTmpl.h>   // base class template
#include <rpg/fts/compressor/Compressor.h>  // base class argument
#include <pscf/cuda/DeviceArray.h>          // base class argument
#include <pscf/cuda/cudaTypes.h>

namespace Pscf {
   extern template
   class AmIteratorTmpl< Rpg::Compressor<1>, DeviceArray<cudaReal> >;
   extern template
   class AmIteratorTmpl< Rpg::Compressor<2>, DeviceArray<cudaReal> >;
   extern template
   class AmIteratorTmpl< Rpg::Compressor<3>, DeviceArray<cudaReal> >;
}
#endif
