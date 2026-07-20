#ifndef RPG_AM_COMPRESSOR_BASE_H
#define RPG_AM_COMPRESSOR_BASE_H

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
#include <rp/fts/compressor/Compressor.h>   // base class argument
#include <pscf/cuda/CudaTp.h>               // base class argument
#include <pscf/cuda/DeviceArray.h>          // base class argument
#include <pscf/cuda/cudaTypes.h>            // base class argument

namespace Pscf {
   extern template
   class AmIteratorTmpl< Rp::Compressor<1, CudaTp<1> >, DeviceArray<cudaReal> >;
   extern template
   class AmIteratorTmpl< Rp::Compressor<2, CudaTp<2> >, DeviceArray<cudaReal> >;
   extern template
   class AmIteratorTmpl< Rp::Compressor<3, CudaTp<3> >, DeviceArray<cudaReal> >;
}
#endif
