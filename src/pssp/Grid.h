#ifndef PSSP_GRID_H
#define PSSP_GRID_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntVec.h"
#include "UnitCell.h"

namespace Pscf { 
namespace Pssp
{ 

   using namespace Util;

   template <int D>
   class Grid {
 
      void setUnitCell(const UnitCell& unitCell)
      void setDimensions(IntVec<D>);

      /**
      * Get vector of grid dimensions.
      */  
      const IntVec<D>& dimensions() const;

      /**
      * Get a grid dimension in a particular direction.
      */  
      int dimension(int i) const;

      /**
      * Return one array index of a grid position.
      */
      int positionId(IntVec<D> position) const;

      /**
      * Return a grid position shifted to convention FFT grid.
      */
      IntVec<D> shiftToFFT(IntVec<D> k) const;

      /**
      * Return a k vector shifted to first Brillouin zone.
      */
      IntVec<D> shiftToBZ(IntVec<D> k) const;
      
   private:
   
      IntVec<D> dimensions_
   
   }

}
}
#endif
