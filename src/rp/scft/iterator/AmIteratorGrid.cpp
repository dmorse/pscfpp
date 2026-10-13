/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions
namespace Pscf {

   // Base class template specializations
   template class
   AmIteratorTmpl< Rp::Iterator<1,CPT>, DeviceArray<double,CPT> >;
   template class
   AmIteratorTmpl< Rp::Iterator<2,CPT>, DeviceArray<double,CPT> >;
   template class
   AmIteratorTmpl< Rp::Iterator<3,CPT>, DeviceArray<double,CPT> >;

   namespace Rp {
      template class AmIteratorGrid<1,CPT>;
      template class AmIteratorGrid<2,CPT>;
      template class AmIteratorGrid<3,CPT>;
   }

}
