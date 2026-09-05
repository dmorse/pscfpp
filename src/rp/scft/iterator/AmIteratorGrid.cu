/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/RField.h>
#include <pscf/backend/cuda/DeviceArray.h>

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>

#include <pscf/backend/cuda/CUT.h>
#include <rp/scft/iterator/AmIteratorGrid.tpp> // template implementation

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
