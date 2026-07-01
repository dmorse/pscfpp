#ifndef RPG_PERTURBATION_FACTORY_H
#define RPG_PERTURBATION_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/perturbation/PerturbationFactory.h>
#include <rpg/system/Types.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PerturbationFactory<1, Rpg::Types<1> >;
      extern template class PerturbationFactory<2, Rpg::Types<2> >;
      extern template class PerturbationFactory<3, Rpg::Types<3> >;
   }
}
#endif
