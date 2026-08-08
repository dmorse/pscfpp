/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>
#include <rp/fts/compressor/CompressorFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
namespace Rp {
   template class CompressorFactory<1,CUT>;
   template class CompressorFactory<2,CUT>;
   template class CompressorFactory<3,CUT>;
}
}
