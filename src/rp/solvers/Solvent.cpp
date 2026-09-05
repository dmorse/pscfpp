/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>

#include <rp/solvers/Solvent.tpp>  // class implementation
#include <pscf/backend/cpp/CPT.h>      // backend identifier class

namespace Pscf {
   namespace Rp {
      template class Solvent<1,CPT>;
      template class Solvent<2,CPT>;
      template class Solvent<3,CPT>;
   }
}
