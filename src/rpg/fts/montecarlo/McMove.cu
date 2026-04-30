/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <rp/fts/montecarlo/McMove.tpp>     // base class implementation

namespace Pscf {
   namespace Rpg {
   
      /*
      * Constructor.
      */
      template <int D>
      McMove<D>::McMove(McSimulator<D>& simulator)
       : Rp::McMove<D, Types<D> >(simulator)
      {}
   
   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMove<1, Rpg::Types<1> >;
      template class McMove<2, Rpg::Types<2> >;
      template class McMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class McMove<1>;
      template class McMove<2>;
      template class McMove<3>;
   }
}
