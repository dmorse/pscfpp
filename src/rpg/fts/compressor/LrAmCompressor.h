#ifndef RPG_LR_AM_COMPRESSOR_H
#define RPG_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrAmCompressor.h>    // base class template
#include <pscf/cuda/CudaTp.h>                    // base class argument

#include <rpg/fts/compressor/AmCompressorBase.h> // indirect base
#include <rpg/fts/compressor/IntraCorrelation.h> // base member
#include <prdc/field/cuda/RField.h>              // base member
#include <prdc/field/cuda/RFieldDft.h>           // base member


// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LrAmCompressor<1, CudaTp<1> >;
      extern template class LrAmCompressor<2, CudaTp<2> >;
      extern template class LrAmCompressor<3, CudaTp<3> >;
   }
}
#endif
