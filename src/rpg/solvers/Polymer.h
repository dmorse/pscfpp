#ifndef RPG_POLYMER_H
#define RPG_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <pscf/cuda/CudaTp.h>      // base class template parameter

namespace Pscf {
   namespace Rp {
      template <int D, class T> class Propagator;
      template <int D, class T> class Block;
   }
}

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class PolymerTmpl< Rp::Block<1, CudaTp<1> >, Rp::Propagator<1, CudaTp<1> > >;
   extern template 
   class PolymerTmpl< Rp::Block<2, CudaTp<2> >, Rp::Propagator<2, CudaTp<2> > >;
   extern template 
   class PolymerTmpl< Rp::Block<3, CudaTp<3> >, Rp::Propagator<3, CudaTp<3> > >;
   namespace Rp {
      extern template class Polymer<1, CudaTp<1> >;
      extern template class Polymer<2, CudaTp<2> >;
      extern template class Polymer<3, CudaTp<3> >;
   }
}
#endif
