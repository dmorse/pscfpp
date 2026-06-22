#ifndef RPG_MIXTURE_MODIFIER_H
#define RPG_MIXTURE_MODIFIER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.h>  // base class template
#include <rpg/solvers/Mixture.h>         // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class MixtureModifier< 1, Rpg::Types<1> >;
      extern template class MixtureModifier< 2, Rpg::Types<2> >;
      extern template class MixtureModifier< 3, Rpg::Types<3> >;
   }
}
#endif
