/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>
#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/scft/iterator/AmIteratorGrid.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template class 
   AmIteratorTmpl< Rp::Iterator<1,CUT>, DeviceArray<cudaReal> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2,CUT>, DeviceArray<cudaReal> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3,CUT>, DeviceArray<cudaReal> >;
   namespace Rp {
      template class AmIteratorGrid<1,CUT>;
      template class AmIteratorGrid<2,CUT>;
      template class AmIteratorGrid<3,CUT>;
   }
}
