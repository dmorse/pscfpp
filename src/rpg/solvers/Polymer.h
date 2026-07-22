#ifndef RPG_POLYMER_H
#define RPG_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <pscf/backends/CUT.h>      // base class template parameter

namespace Pscf {
   namespace Rp {
      template <int D, class T> class Propagator;
      template <int D, class T> class Block;
   }
}

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class PolymerTmpl< Rp::Block<1,CUT>, Rp::Propagator<1,CUT> >;
   extern template 
   class PolymerTmpl< Rp::Block<2,CUT>, Rp::Propagator<2,CUT> >;
   extern template 
   class PolymerTmpl< Rp::Block<3,CUT>, Rp::Propagator<3,CUT> >;
   namespace Rp {
      extern template class Polymer<1,CUT>;
      extern template class Polymer<2,CUT>;
      extern template class Polymer<3,CUT>;
   }
}
#endif
