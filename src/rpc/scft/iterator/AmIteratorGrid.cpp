/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"                    // class header

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <rpc/field/Mask.h>

#include <pscf/cpu/FftwDRArray.h>
#include <prdc/field/cpu/RField.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/iterator/AmIteratorGrid.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {

   template class 
   AmIteratorTmpl< Rp::Iterator<1, CppTp<1> >, FftwDRArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2, CppTp<2> >, FftwDRArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3, CppTp<3> >, FftwDRArray<double> >;

   namespace Rp {
      template class AmIteratorGrid<1, CppTp<1> >;
      template class AmIteratorGrid<2, CppTp<2> >;
      template class AmIteratorGrid<3, CppTp<3> >;
   }

}
