/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"              // class header
#include <rpg/field/FieldIo.tpp> 
#include <prdc/field/cuda/RField.tpp> 

#include <rp/field/CFields.tpp>   // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class CFields<1, CudaTp<1> >;
      template class CFields<2, CudaTp<2> >;
      template class CFields<3, CudaTp<3> >;
   } 
}
