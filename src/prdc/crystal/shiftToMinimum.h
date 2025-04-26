#ifndef PRDC_SHIFT_MINIMUM_H
#define PRDC_SHIFT_MINIMUM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <iostream>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Returns minimum magnitude image of DFT wavevector.
   *
   * \param v IntVec<D> containing integer indices of wavevector.
   * \param d dimensions of the discrete Fourier transform grid.
   * \param cell UnitCell
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   IntVec<D> shiftToMinimum(IntVec<D>& v, IntVec<D> d, UnitCell<D> const & cell);

   // Explicit specializations
   template <> 
   IntVec<1> shiftToMinimum(IntVec<1>& v, IntVec<1> d, UnitCell<1> const & cell);

   template <> 
   IntVec<2> shiftToMinimum(IntVec<2>& v, IntVec<2> d, UnitCell<2> const & cell);

   template <> 
   IntVec<3> shiftToMinimum(IntVec<3>& v, IntVec<3> d, UnitCell<3> const & cell);

}
}
#endif
