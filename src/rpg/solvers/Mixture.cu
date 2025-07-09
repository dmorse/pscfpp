/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"

namespace Pscf {
   namespace Prdc { 
      template class MixtureReal<1, Rpg::Polymer<1>, Rpg::Solvent<1> >;
      template class MixtureReal<2, Rpg::Polymer<2>, Rpg::Solvent<2> >;
      template class MixtureReal<3, Rpg::Polymer<3>, Rpg::Solvent<3> >;
   }
   namespace Rpg { 
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}
