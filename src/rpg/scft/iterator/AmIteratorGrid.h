#ifndef RPG_AM_ITERATOR_GRID_H
#define RPG_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.h>  // class template
#include <pscf/backends/CUT.h>                 // class argument
#include <rpg/scft/iterator/Iterator.h>       // indirect base class 
#include <pscf/cuda/DeviceArray.h>            // base class argument

// Explicit instantiation declarations
namespace Pscf {

   extern template class 
   AmIteratorTmpl<Rp::Iterator<1,CUT>, DeviceArray<cudaReal> >;
   extern template class 
   AmIteratorTmpl<Rp::Iterator<2,CUT>, DeviceArray<cudaReal> >;
   extern template class 
   AmIteratorTmpl<Rp::Iterator<3,CUT>, DeviceArray<cudaReal> >;

   namespace Rp {
      extern template class AmIteratorGrid<1,CUT>;
      extern template class AmIteratorGrid<2,CUT>;
      extern template class AmIteratorGrid<3,CUT>;
   }

} 
#endif
