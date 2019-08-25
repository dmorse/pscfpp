#ifndef PSSP_FIELD_UTIL_H
#define PSSP_FIELD_UTIL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp/basis/Basis.h>
#include <pssp/field/RFieldDft.h>
#include <util/containers/DArray.h>

namespace Pscf { 
namespace Pssp 
{ 

   using namespace Util;

   /**
   * Convert field from symmetry-adapted Fourier transform to DFT (k-grid).
   *
   * \param basis symmetry adapted discrete Fourier basis
   * \param components coefficients of symmetry-adapted basis functions
   * \param dft discrete Fourier transform of a real field
   * \ingroup Pssp_Field_Module
   */
   template <int D>
   void convertFieldBasisToDft(Basis<D> const & basis,
                               DArray<double> const& components, 
                               RFieldDft<D>& dft);

   /**
   * Convert DFT of field (k-grid) to symmetry-adapted Fourier transform.
   *
   * \param basis symmetry adapted discrete Fourier basis
   * \param dft complex DFT representation of a field.
   * \param components coefficients of symmetry-adapted basis functions.
   * \ingroup Pssp_Field_Module
   */
   template <int D>
   void convertFieldDftToBasis(Basis<D> const & basis,
                               RFieldDft<D>& dft, 
                               DArray<double>& components);

} // namespace Pscf:Pssp
} // namespace Pscf
#include "fieldUtil.tpp"
#endif
