/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"                    // header
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <rp/fts/analyzer/HamiltonianAnalyzer.tpp>  // implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(Simulator<D>& simulator,
                                   System<D>& system)
    : Rp::HamiltonianAnalyzer< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class HamiltonianAnalyzer< 1, Rpg::Types<1> >;
      template class HamiltonianAnalyzer< 2, Rpg::Types<2> >;
      template class HamiltonianAnalyzer< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class HamiltonianAnalyzer<1>;
      template class HamiltonianAnalyzer<2>;
      template class HamiltonianAnalyzer<3>;
   }
}
