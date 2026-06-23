/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>

#include <rp/fts/analyzer/AverageAnalyzer.tpp>

namespace Pscf {
   namespace Rpc {

      /// Constructor.
      template <int D>
      AverageAnalyzer<D>::AverageAnalyzer(Simulator<D>& simulator,
		                          Rp::System<D, Rpc::Types<D> >& system)
       : Rp::AverageAnalyzer< D, Types<D> >(simulator, system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class
      AverageAnalyzer<1, Rpc::Types<1> >;
      template class
      AverageAnalyzer<2, Rpc::Types<2> >;
      template class
      AverageAnalyzer<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class AverageAnalyzer<1>;
      template class AverageAnalyzer<2>;
      template class AverageAnalyzer<3>;
   }
}
