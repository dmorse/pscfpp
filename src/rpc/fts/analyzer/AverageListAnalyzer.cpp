/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"                    // header
#include <rpc/system/System.h>
#include <rp/fts/analyzer/AverageListAnalyzer.tpp>  // implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   AverageListAnalyzer<D>::AverageListAnalyzer(Rp::Simulator<D, Rpc::Types<D> >& simulator,
                                   Rp::System<D, Rpc::Types<D> >& system)
    : Rp::AverageListAnalyzer< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageListAnalyzer<1, Rpc::Types<1> >;
      template class AverageListAnalyzer<2, Rpc::Types<2> >;
      template class AverageListAnalyzer<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class AverageListAnalyzer<1>;
      template class AverageListAnalyzer<2>;
      template class AverageListAnalyzer<3>;
   }
}
