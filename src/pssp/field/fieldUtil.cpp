/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "fieldUtil.tpp"

namespace Pscf {
namespace Pssp {

   template <>
   void convertFieldBasisToDft(Basis<3> const & basis,
                               DArray<double> const& components, 
                               RFieldDft<3>& dft);

   template <>
   void convertFieldDftToBasis(Basis<3> const & basis,
                               RFieldDft<3> const & dft, 
                               DArray<double>& components);

} // namespace Pscf::Pssp
} // namespace Pscf
