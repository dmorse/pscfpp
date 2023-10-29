/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/Field.tpp>

namespace Pscf { 
namespace Prdc { 
namespace Cuda { 

   template class Field<cudaReal>;
   template class Field<cudaComplex>;

}
}
}
