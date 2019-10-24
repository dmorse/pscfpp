#ifndef PSPG_SHIFT_MINIMUM_H
#define PSPG_SHIFT_MINIMUM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <iostream>

namespace Pscf{
namespace Pspg{

   using namespace Util;

   /**
   * Returns minimum magnitude image of DFT wavevector.
   *
   * \param v IntVec<D> containing integer indices of wavevector.
   * \param d dimensions of the discrete Fourier transform grid.
   * \param cell UnitCell
   */
   template <int D>
   IntVec<D> shiftToMinimum(IntVec<D>& v, const IntVec<D> d, const UnitCell<D> cell);

   // Explicit specializations
   // The explicit specializations assumes that the value of IntVec is strictly
   // non-negative. This is always safe if IntVec<D>& v is provided 
   // by MeshIterator
   template <> 
   IntVec<1> shiftToMinimum(IntVec<1>& v, const IntVec<1> d, const UnitCell<1> cell);

   template <> 
   IntVec<2> shiftToMinimum(IntVec<2>& v, const IntVec<2> d, const UnitCell<2> cell);

   template <> 
   IntVec<3> shiftToMinimum(IntVec<3>& v, const IntVec<3> d, const UnitCell<3> cell);


}
}
#endif
