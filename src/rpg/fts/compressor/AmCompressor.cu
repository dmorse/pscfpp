/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/fts/compressor/AmCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template 
      class AmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      template 
      class AmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      template 
      class AmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
}
