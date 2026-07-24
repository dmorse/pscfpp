/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>  // backend type class

//#include <prdc/field/cpu/RField.h>
//#include <pscf/mesh/Mesh.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/solvers/Solvent.tpp>

namespace Pscf {
   namespace Rp {
      template class Solvent<1,CPT>;
      template class Solvent<2,CPT>;
      template class Solvent<3,CPT>;
   }
}
