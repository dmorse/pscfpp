/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>

#include <rp/fts/analyzer/AverageAnalyzer.tpp>

namespace Pscf {
   namespace Rpg {

      /// Constructor.
      template <int D>
      AverageAnalyzer<D>::AverageAnalyzer(Rp::Simulator<D, Rpg::Types<D> >& simulator,
		                          Rp::System<D, Rpg::Types<D> >& system)
       : Rp::AverageAnalyzer< D, Types<D> >(simulator, system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class
      AverageAnalyzer<1, Rpg::Types<1> >;
      template class
      AverageAnalyzer<2, Rpg::Types<2> >;
      template class
      AverageAnalyzer<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class AverageAnalyzer<1>;
      template class AverageAnalyzer<2>;
      template class AverageAnalyzer<3>;
   }
}
