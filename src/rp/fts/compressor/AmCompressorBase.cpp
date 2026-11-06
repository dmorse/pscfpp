/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressorBase.h"

#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template 
   class AmIteratorTmpl< Rp::Compressor<1,CPT>, BackendArray<double,CPT> >;
   template 
   class AmIteratorTmpl< Rp::Compressor<2,CPT>, BackendArray<double,CPT> >; 
   template 
   class AmIteratorTmpl< Rp::Compressor<3,CPT>, BackendArray<double,CPT> >;
}
