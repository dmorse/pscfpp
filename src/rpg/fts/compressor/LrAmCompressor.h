#ifndef RPG_LR_AM_COMPRESSOR_H
#define RPG_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrAmCompressor.h>    // base class template
#include <rpg/system/Types.h>                    // base class argument
#include <rpg/fts/compressor/IntraCorrelation.h> // base member
#include <prdc/field/cuda/RField.h>              // base member
#include <prdc/field/cuda/RFieldDft.h>           // base member
#include <rpg/fts/compressor/AmCompressorBase.h> // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template
      class LrAmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      extern template
      class LrAmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      extern template
      class LrAmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
}
#endif
