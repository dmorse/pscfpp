/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"              // class header
#include <rpg/field/FieldIo.tpp> 
#include <prdc/cuda/RField.tpp> 

#include <rp/field/CFields.tpp>   // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      template class CFields<1, Rpg::Types<1> >;
      template class CFields<2, Rpg::Types<2> >;
      template class CFields<3, Rpg::Types<3> >;
   } 
   namespace Rpg {
      template class CFields<1>;
      template class CFields<2>;
      template class CFields<3>;
   } 
}
