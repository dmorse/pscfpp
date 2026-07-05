#ifndef RPG_LR_COMPRESSOR_H
#define RPG_LR_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/LrCompressor.h>      // direct base template
#include <rpg/system/Types.h>                    // direct base argument
#include <rpg/fts/compressor/IntraCorrelation.h> // direct base member
#include <prdc/field/cuda/RField.h>                    // direct base member
#include <prdc/field/cuda/RFieldDft.h>                 // direct base member
#include <rp/fts/compressor/Compressor.h>        // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LrCompressor<1, Rpg::Types<1> >;
      extern template class LrCompressor<2, Rpg::Types<2> >;
      extern template class LrCompressor<3, Rpg::Types<3> >;
   }
}
#endif
