/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* Instantiated specializations of AmIteratorTmpl used as base classes 
* for specializations of Rpg::AmCompressor and Rpg::LrAmCompressor.
*/

#include "AmCompressorBase.h"
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template 
   class AmIteratorTmpl< Rp::Compressor<1, Rpg::Types<1> >, DeviceArray<cudaReal> >;
   template 
   class AmIteratorTmpl< Rp::Compressor<2, Rpg::Types<2> >, DeviceArray<cudaReal> >; 
   template 
   class AmIteratorTmpl< Rp::Compressor<3, Rpg::Types<3> >, DeviceArray<cudaReal> >;
}
