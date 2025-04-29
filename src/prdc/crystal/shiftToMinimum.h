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
   * \ingroup Prdc_Crystal_Module
   *
   * \param v  IntVec<D> containing integer indices of a wavevector
   * \param d  dimensions of the discrete Fourier transform grid
   * \param cell  UnitCell
   * \return  integer indices of minimum image of v
   */
   template <int D>
   IntVec<D> shiftToMinimum(IntVec<D> const & v, 
                            IntVec<D> const & d, 
                            UnitCell<D> const & cell);

   // Explicit specializations for D=1, 2 and 3

   template <> 
   IntVec<1> shiftToMinimum(IntVec<1> const & v, 
                            IntVec<1> const & d, 
                            UnitCell<1> const & cell);

   template <> 
   IntVec<2> shiftToMinimum(IntVec<2> const & v, 
                            IntVec<2> const & d, 
                            UnitCell<2> const & cell);

   template <> 
   IntVec<3> shiftToMinimum(IntVec<3> const & v, 
                            IntVec<3> const & d, 
                            UnitCell<3> const & cell);

}
}
#endif
