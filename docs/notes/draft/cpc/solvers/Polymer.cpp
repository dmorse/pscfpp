/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.tpp"

namespace Pscf {

   template 
   class PolymerTmpl< Cpc::Block<1>, Cpc::Propagator<1>, std::complex<double> >;
   template 
   class PolymerTmpl< Cpc::Block<2>, Cpc::Propagator<2>, std::complex<double> >;
   template 
   class PolymerTmpl< Cpc::Block<3>, Cpc::Propagator<3>, std::complex<double> >;

   namespace Cpc {
      template class Polymer<1>;
      template class Polymer<2>;
      template class Polymer<3>;
   }

}
