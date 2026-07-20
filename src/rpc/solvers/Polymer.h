#ifndef RPC_POLYMER_H
#define RPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <pscf/cpu/CppTp.h>      // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class PolymerTmpl< Rp::Block<1, CppTp<1> >, Rp::Propagator<1, CppTp<1> > >;
   extern template 
   class PolymerTmpl< Rp::Block<2, CppTp<2> >, Rp::Propagator<2, CppTp<2> > >;
   extern template 
   class PolymerTmpl< Rp::Block<3, CppTp<3> >, Rp::Propagator<3, CppTp<3> > >;
   namespace Rp {
      extern template class Polymer<1, CppTp<1> >;
      extern template class Polymer<2, CppTp<2> >;
      extern template class Polymer<3, CppTp<3> >;
   }
}
#endif
