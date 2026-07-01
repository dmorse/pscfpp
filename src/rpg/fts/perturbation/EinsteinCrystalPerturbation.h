#ifndef RPG_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RPG_EINSTEIN_CRYSTAL_PERTURBATION_H

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.h> // base class
#include <rpg/system/Types.h>                               // base argument
#include <rpg/fts/perturbation/Perturbation.h>              // indirect base
#include <prdc/field/cuda/RField.h>                         // base member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class EinsteinCrystalPerturbation<1, Rpg::Types<1> >;
      extern template class EinsteinCrystalPerturbation<2, Rpg::Types<2> >;
      extern template class EinsteinCrystalPerturbation<3, Rpg::Types<3> >;
   }
}
#endif
