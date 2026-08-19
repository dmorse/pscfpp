/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <pscf/backends/CPT.h>     // backend type class
#include <rp/solvers/Solvent.tpp>  // class implementation

namespace Pscf {
   namespace Rp {
      template class Solvent<1,CPT>;
      template class Solvent<2,CPT>;
      template class Solvent<3,CPT>;
   }
}
