/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

/*
* Instantiated specializations of AmIteratorTmpl used as base classes 
* for specializations of Rpg::AmCompressor and Rpg::LrAmCompressor.
*/

#include "AmCompressorBase.h"

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template 
   class AmIteratorTmpl< Rp::Compressor<1,CUT>, BackendArray<cudaReal,CUT> >;
   template 
   class AmIteratorTmpl< Rp::Compressor<2,CUT>, BackendArray<cudaReal,CUT> >; 
   template 
   class AmIteratorTmpl< Rp::Compressor<3,CUT>, BackendArray<cudaReal,CUT> >;
}
