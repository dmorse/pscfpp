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
namespace Rpg {

   // Constructor.
   template <int D>
   AnalyzerManager<D>::AnalyzerManager(Rp::Simulator<D, Rpg::Types<D> >& simulator,
                                       Rp::System<D, Rpg::Types<D> >& system)
    : Rp::AnalyzerManager<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class AnalyzerManager<1, Rpg::Types<1> >;
      template class AnalyzerManager<2, Rpg::Types<2> >;
      template class AnalyzerManager<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class AnalyzerManager<1>;
      template class AnalyzerManager<2>;
      template class AnalyzerManager<3>;
   }
}
