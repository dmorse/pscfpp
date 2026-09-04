/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>
#include <pscf/backend/cpp/CpuVecRandom.h>
#include <pscf/backend/CPT.h>

#include <rp/fts/simulator/Simulator.tpp>  // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Simulator<1,CPT>;
      template class Simulator<2,CPT>;
      template class Simulator<3,CPT>;
   }
}
