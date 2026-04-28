/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"                     // class header
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <rpc/field/Mask.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/ScftThermo.tpp>           // base class implementation

namespace Pscf {
   namespace Rpc {

      /*
      * Constructor.
      */
      template <int D>
      ScftThermo<D>::ScftThermo(System<D> const & system)
       : Rp::ScftThermo<D, Types<D> >(system)
      {};

   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1, Rpc::Types<1> >;
      template class ScftThermo<2, Rpc::Types<2> >;
      template class ScftThermo<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;
   }
}
