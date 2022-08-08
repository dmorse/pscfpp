/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIterator.tpp"
#include "AmIterator.tpp"

namespace Pscf {
namespace Pspc {

   template class FilmIterator<1, AmIterator<1> >;
   template class FilmIterator<2, AmIterator<2> >;
   template class FilmIterator<3, AmIterator<3> >;

}
} 
