/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldState.tpp"

namespace Pscf {
   namespace Rpc
   {

      template class FieldState< 1, DArray<double> >;
      template class FieldState< 2, DArray<double> >;
      template class FieldState< 3, DArray<double> >;

   } // namespace Rpc
} // namespace Pscf
