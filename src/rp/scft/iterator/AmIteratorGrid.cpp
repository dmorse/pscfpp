/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>
#include <prdc/field/RField.h>

#include <rp/scft/iterator/AmIteratorGrid.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {

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
