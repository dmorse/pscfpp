/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/system/System.h>
#include "BasisFieldState.tpp"

// Explicit instantiation declarations
namespace Pscf {
   namespace Rpc {
      template class FieldState< 1, DArray<double>, System<1> >;
      template class FieldState< 2, DArray<double>, System<2> >;
      template class FieldState< 3, DArray<double>, System<3> >;
      template class BasisFieldState<1>;
      template class BasisFieldState<2>;
      template class BasisFieldState<3>;
   }
}
