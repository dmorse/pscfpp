/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EinsteinCrystalPerturbation.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.tpp>

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(
                                                  Simulator<D>& simulator)
    : Rp::EinsteinCrystalPerturbation<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EinsteinCrystalPerturbation<1, Rpc::Types<1> >;
      template class EinsteinCrystalPerturbation<2, Rpc::Types<2> >;
      template class EinsteinCrystalPerturbation<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class EinsteinCrystalPerturbation<1>;
      template class EinsteinCrystalPerturbation<2>;
      template class EinsteinCrystalPerturbation<3>;
   }
}
