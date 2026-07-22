#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_H

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.h> // base class
#include <pscf/backends/CPT.h>                                // base argument
#include <rpc/fts/perturbation/Perturbation.h>               // indirect base
#include <prdc/field/cpu/RField.h>                           // base member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class EinsteinCrystalPerturbation<1,CPT>;
      extern template 
      class EinsteinCrystalPerturbation<2,CPT>;
      extern template 
      class EinsteinCrystalPerturbation<3,CPT>;
   }
}
#endif
