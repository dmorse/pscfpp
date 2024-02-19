/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIterator.tpp"
#include "AmIteratorBasis.tpp"

namespace Pscf {
namespace Rpc {

   template class FilmIterator<1, AmIteratorBasis<1> >;
   template class FilmIterator<2, AmIteratorBasis<2> >;
   template class FilmIterator<3, AmIteratorBasis<3> >;

}
} 
