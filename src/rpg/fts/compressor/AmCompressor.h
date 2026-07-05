#ifndef RPG_AM_COMPRESSOR_H
#define RPG_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/AmCompressor.h>       // direct base template
#include <rpg/system/Types.h>                     // direct base argument
#include <rpg/fts/compressor/AmCompressorBase.h>  // indirect base 
#include <prdc/field/cuda/RField.h>               // direct base member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class AmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      extern template 
      class AmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      extern template 
      class AmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
}
#endif
