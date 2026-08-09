#ifndef RPC_LR_COMPRESSOR_H
#define RPC_LR_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrCompressor.h>       // base class template
#include <pscf/backends/CPT.h>                     // base class argument
#include <rp/fts/compressor/IntraCorrelation.h>  // base member
#include <prdc/field/cpu/RField.h>                // base member
#include <prdc/field/cpu/RFieldDft.h>             // base member
#include <rp/fts/compressor/Compressor.h>         // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LrCompressor<1,CPT>;
      extern template class LrCompressor<2,CPT>;
      extern template class LrCompressor<3,CPT>;
   }
}
#endif
