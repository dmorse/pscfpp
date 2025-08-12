/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.tpp"

namespace Pscf {
   namespace Rpc
   {

      template class SweepFactory<1>;
      template class SweepFactory<2>;
      template class SweepFactory<3>;

   } // namespace Rpc
} // namespace Pscf
