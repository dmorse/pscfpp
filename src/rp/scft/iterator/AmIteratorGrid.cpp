/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/RField.h>
#include <pscf/cpu/FftwDRArray.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <pscf/backends/CPT.h>
#include <rp/scft/iterator/AmIteratorGrid.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {

   template class 
   AmIteratorTmpl< Rp::Iterator<1,CPT>, FftwDRArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2,CPT>, FftwDRArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3,CPT>, FftwDRArray<double> >;

   namespace Rp {
      template class AmIteratorGrid<1,CPT>;
      template class AmIteratorGrid<2,CPT>;
      template class AmIteratorGrid<3,CPT>;
   }

}
