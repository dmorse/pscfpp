/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.tpp"

namespace Pscf {

   template class PolymerTmpl< Rpc::Block<1>, Rpc::Propagator<1> >;
   template class PolymerTmpl< Rpc::Block<2>, Rpc::Propagator<2> >;
   template class PolymerTmpl< Rpc::Block<3>, Rpc::Propagator<3> >;

   namespace Rpc {
      template class Polymer<1>;
      template class Polymer<2>;
      template class Polymer<3>;
   }

}
