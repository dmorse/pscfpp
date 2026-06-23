/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/Analyzer.tpp>

namespace Pscf {
   namespace Rpc {

      /// Constructor.
      template <int D>
      Analyzer<D>::Analyzer(Simulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system)
       : Rp::Analyzer<D, Simulator<D>, Rp::System<D, Rpc::Types<D> > >(simulator, system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Analyzer<1, Rpc::Simulator<1>, Rp::System<1, Rpc::Types<1> > >;
      template class Analyzer<2, Rpc::Simulator<2>, Rp::System<2, Rpc::Types<2> > >;
      template class Analyzer<3, Rpc::Simulator<3>, Rp::System<3, Rpc::Types<3> > >;
   } 
   namespace Rpc {
      template class Analyzer<1>;
      template class Analyzer<2>;
      template class Analyzer<3>;
   } 
} 
