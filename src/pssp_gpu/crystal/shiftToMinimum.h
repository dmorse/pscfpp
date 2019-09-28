#ifndef PSSP_GPU_SHIFT_MINIMUM_H
#define PSSP_GPU_SHIFT_MINIMUM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>
#include <iostream>

namespace Pscf{
namespace Pssp_gpu{

   using namespace Util;

   /**
   * Returns minimum magnitude image of DFT wavevector.
   *
   * \param v IntVec<D> containing integer indices of wavevector.
   * \param d dimensions of the discrete Fourier transform grid.
   * \param cell UnitCell
   */
   template <int D>
   IntVec<D> shiftToMinimum(IntVec<D>& v, const IntVec<D> d);

   template <int D>
   void shiftToMinimum(const IntVec<D>& v, const IntVec<D>& d, int* waveBz);
   // Explicit specializations
   // The explicit specializations assumes that the value of IntVec is strictly
   // non-negative. This is always safe if IntVec<D>& v is provided 
   // by MeshIterator

   template <int D> 
   IntVec<D> shiftToMinimum(IntVec<D>& v, const IntVec<D> d)
   {
      IntVec<D> u;
      for( int i = 0; i < D; ++i) {
         if (v[i] > d[i]/2) {
            u[i] = v[i] - d[i];
         } else {
            u[i] = v[i];
         }
      }
      return u;
   }

   template<int D>
   void shiftToMinimum(const IntVec<D>& v, const IntVec<D>& d, int* waveBz) {      
      for(int i = 0; i < D; i++) {
         if (v[i] > d[i] / 2 ) {
            waveBz[i] = v[i] - d[i];
         } else {
            waveBz[i] = v[i];
         }
      }
   }

}
}
#endif
