/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"

// Explicit instantiation definitions
namespace Pscf {
   namespace Cp { 
      template 
      class Mixture<1, Cpc::Polymer<1>, Cpc::Solvent<1>, Cpc::Types<1> >;
      template 
      class Mixture<2, Cpc::Polymer<2>, Cpc::Solvent<2>, Cpc::Types<2> >;
      template 
      class Mixture<3, Cpc::Polymer<3>, Cpc::Solvent<3>, Cpc::Types<3> >;
   }
   namespace Cpc { 
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}
