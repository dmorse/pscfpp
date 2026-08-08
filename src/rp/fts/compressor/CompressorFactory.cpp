/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>
#include <rp/fts/compressor/CompressorFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
namespace Rp {
   template class CompressorFactory<1,CPT>;
   template class CompressorFactory<2,CPT>;
   template class CompressorFactory<3,CPT>;
}
}
