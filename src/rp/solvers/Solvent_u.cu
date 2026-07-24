/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>    // backend type class

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/solvers/Solvent.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Solvent<1,CUT>;
      template class Solvent<2,CUT>;
      template class Solvent<3,CUT>;
   }
}
