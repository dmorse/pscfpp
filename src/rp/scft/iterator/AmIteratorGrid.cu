/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation definitions
namespace Pscf {

   template class
   AmIteratorTmpl< Rp::Iterator<1,CUT>, DeviceArray<cudaReal,CUT> >;
   template class
   AmIteratorTmpl< Rp::Iterator<2,CUT>, DeviceArray<cudaReal,CUT> >;
   template class
   AmIteratorTmpl< Rp::Iterator<3,CUT>, DeviceArray<cudaReal,CUT> >;

   namespace Rp {
      template class AmIteratorGrid<1,CUT>;
      template class AmIteratorGrid<2,CUT>;
      template class AmIteratorGrid<3,CUT>;
   }

}
