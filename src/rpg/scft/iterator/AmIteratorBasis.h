#ifndef RPG_AM_ITERATOR_BASIS_H
#define RPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.h>  // direct base class 
#include <pscf/backends/CUT.h>                  // direct base argument
#include <rpg/scft/iterator/Iterator.h>        // indirect base argument
#include <util/containers/DArray.h>            // indirect base argument

// Explicit instantiation declarations
namespace Pscf {
   extern template class AmIteratorTmpl<Rp::Iterator<1,CUT>, DArray<double> >;
   extern template class AmIteratorTmpl<Rp::Iterator<2,CUT>, DArray<double> >;
   extern template class AmIteratorTmpl<Rp::Iterator<3,CUT>, DArray<double> >;
   namespace Rp {
      extern template class AmIteratorBasis<1,CUT>;
      extern template class AmIteratorBasis<2,CUT>;
      extern template class AmIteratorBasis<3,CUT>;
   } 
}
#endif
