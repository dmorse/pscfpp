/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/system/System.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <rp/fts/montecarlo/McMove.tpp>     // base class implementation

namespace Pscf {
   namespace Rpc {
   
      /*
      * Constructor.
      */
      template <int D>
      McMove<D>::McMove(Rp::McSimulator<D, Rpc::Types<D> >& simulator)
       : Rp::McMove<D, Types<D> >(simulator)
      {}
   
   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMove<1, Rpc::Types<1> >;
      template class McMove<2, Rpc::Types<2> >;
      template class McMove<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class McMove<1>;
      template class McMove<2>;
      template class McMove<3>;
   }
}
