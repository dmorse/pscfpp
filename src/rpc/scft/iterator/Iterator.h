#ifndef RPC_ITERATOR_H
#define RPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/Iterator.h>  // base class template
#include <pscf/cpu/CppTp.h>           // base class template

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Iterator<1, CppTp<1> >;
      extern template class Iterator<2, CppTp<2> >;
      extern template class Iterator<3, CppTp<3> >;
   }
} 
#endif
