/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerManager.h"                   // class header
#include "AnalyzerFactory.h"
#include <rp/fts/analyzer/AnalyzerManager.tpp> // base class implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   AnalyzerManager<D>::AnalyzerManager(Rp::Simulator<D, Rpc::Types<D> >& simulator,
                                       Rp::System<D, Rpc::Types<D> >& system)
    : Rp::AnalyzerManager<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class AnalyzerManager<1, Rpc::Types<1> >;
      template class AnalyzerManager<2, Rpc::Types<2> >;
      template class AnalyzerManager<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class AnalyzerManager<1>;
      template class AnalyzerManager<2>;
      template class AnalyzerManager<3>;
   }
}
