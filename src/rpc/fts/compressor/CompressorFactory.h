#ifndef RPC_COMPRESSOR_FACTORY_H
#define RPC_COMPRESSOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/CompressorFactory.h>  // base class template
#include <pscf/backends/CPT.h>                     // base class argument

// Explicit instantiation declarations
namespace Pscf {
namespace Rp {
   extern template class CompressorFactory<1,CPT>;
   extern template class CompressorFactory<2,CPT>;
   extern template class CompressorFactory<3,CPT>;
}
}
#endif
