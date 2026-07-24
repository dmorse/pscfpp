/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/Compressor.tpp>
#include <pscf/backends/CUT.h>

// Explicit instantiation definitions
namespace Pscf {
namespace Rp {
   template class Compressor<1,CUT>;
   template class Compressor<2,CUT>;
   template class Compressor<3,CUT>;
}
}
