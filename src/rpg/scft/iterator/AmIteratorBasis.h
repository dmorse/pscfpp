#ifndef RPG_AM_ITERATOR_BASIS_H
#define RPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.h>  // direct base class 
#include <rpg/system/Types.h>                  // direct base argument
#include <rpg/scft/iterator/Iterator.h>        // indirect base argument
#include <util/containers/DArray.h>            // indirect base argument

// Explicit instantiation declarations
namespace Pscf {
   extern template class AmIteratorTmpl<Rp::Iterator<1, Rpg::Types<1> >, DArray<double> >;
   extern template class AmIteratorTmpl<Rp::Iterator<2, Rpg::Types<2> >, DArray<double> >;
   extern template class AmIteratorTmpl<Rp::Iterator<3, Rpg::Types<3> >, DArray<double> >;
   namespace Rp {
      extern template class AmIteratorBasis<1, Rpg::Types<1> >;
      extern template class AmIteratorBasis<2, Rpg::Types<2> >;
      extern template class AmIteratorBasis<3, Rpg::Types<3> >;
   } 
}
#endif
