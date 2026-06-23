/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Perturbation.h"
#include <rpc/fts/simulator/Simulator.h>
#include <prdc/cpu/RField.h>

#include <rp/fts/perturbation/Perturbation.tpp> // base class implementation

namespace Pscf {
namespace Rpc {

   /*
   * Constructor.
   */
   template <int D>
   Perturbation<D>::Perturbation(Rp::Simulator<D, Rpc::Types<D> >& simulator)
    : Rp::Perturbation<D, Types<D> > (simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Perturbation<1, Rpc::Types<1> >;
      template class Perturbation<2, Rpc::Types<2> >;
      template class Perturbation<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class Perturbation<1>;
      template class Perturbation<2>;
      template class Perturbation<3>;
   }
}
