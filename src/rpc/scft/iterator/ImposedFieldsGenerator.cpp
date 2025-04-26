/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ImposedFieldsGenerator.tpp"

namespace Pscf {
namespace Rpc {
   template class ImposedFieldsGenerator<1>;
   template class ImposedFieldsGenerator<2>;
   template class ImposedFieldsGenerator<3>;
}
}