/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EinsteinCrystalPerturbation.h"
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <rpg/field/WFields.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.tpp>

#if 0
namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(
                                                  Rp::Simulator<D, Rpg::Types<D> >& simulator)
    : Rp::EinsteinCrystalPerturbation<D, Types<D> >(simulator)
   {}

}
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EinsteinCrystalPerturbation<1, Rpg::Types<1> >;
      template class EinsteinCrystalPerturbation<2, Rpg::Types<2> >;
      template class EinsteinCrystalPerturbation<3, Rpg::Types<3> >;
   }
}
