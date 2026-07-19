/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenExt.tpp"

namespace Pscf {
namespace Rp {

   template class FilmFieldGenExt<1, Rpc::Types<1> >;
   template class FilmFieldGenExt<2, Rpc::Types<2> >;
   template class FilmFieldGenExt<3, Rpc::Types<3> >;

}
}
