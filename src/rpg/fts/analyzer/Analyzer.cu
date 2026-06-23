/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/Analyzer.tpp>

// Subclass constructor definition
namespace Pscf {
   namespace Rpg {

      /*
      * Constructor.
      */
      template <int D>
      Analyzer<D>::Analyzer(Rp::Simulator<D, Rpg::Types<D> >& simulator, Rp::System<D, Rpg::Types<D> >& system)
       : Rp::Analyzer<D, Rp::Simulator<D, Rpg::Types<D> >, Rp::System<D, Rpg::Types<D> > >(simulator, system)
      {}

   }
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class Analyzer<1, Rp::Simulator<1, Rpg::Types<1> >, Rp::System<1, Rpg::Types<1> > >;
      template class Analyzer<2, Rp::Simulator<2, Rpg::Types<2> >, Rp::System<2, Rpg::Types<2> > >;
      template class Analyzer<3, Rp::Simulator<3, Rpg::Types<3> >, Rp::System<3, Rpg::Types<3> > >;
   } 
   namespace Rpg {
      template class Analyzer<1>;
      template class Analyzer<2>;
      template class Analyzer<3>;
   } 
} 
