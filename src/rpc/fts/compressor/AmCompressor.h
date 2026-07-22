#ifndef RPC_AM_COMPRESSOR_H
#define RPC_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/AmCompressor.h>       // base class template
#include <pscf/backends/CPT.h>                     // base class argument
#include <rpc/fts/compressor/AmCompressorBase.h>  // indirect base
#include <prdc/field/cpu/RField.h>                // base member
#include <pscf/cpu/FftwDRArray.h>           // base member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AmCompressor<1,CPT>;
      extern template class AmCompressor<2,CPT>;
      extern template class AmCompressor<3,CPT>;
   }
}
#endif
