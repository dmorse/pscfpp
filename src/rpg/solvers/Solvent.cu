/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/mesh/Mesh.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/solvers/Solvent.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Solvent<1, Rpg::Types<1> >;
      template class Solvent<2, Rpg::Types<2> >;
      template class Solvent<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Solvent<1>;
      template class Solvent<2>;
      template class Solvent<3>;
   }
}
