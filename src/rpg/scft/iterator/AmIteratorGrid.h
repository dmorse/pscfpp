#ifndef RPG_AM_ITERATOR_GRID_H
#define RPG_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.h>  // class template
#include <pscf/cuda/CudaTp.h>                 // class argument
#include <rpg/scft/iterator/Iterator.h>       // indirect base class 
#include <pscf/cuda/DeviceArray.h>            // base class argument

// Explicit instantiation declarations
namespace Pscf {

   extern template class 
   AmIteratorTmpl<Rp::Iterator<1, CudaTp<1> >, DeviceArray<cudaReal> >;
   extern template class 
   AmIteratorTmpl<Rp::Iterator<2, CudaTp<2> >, DeviceArray<cudaReal> >;
   extern template class 
   AmIteratorTmpl<Rp::Iterator<3, CudaTp<3> >, DeviceArray<cudaReal> >;

   namespace Rp {
      extern template class AmIteratorGrid<1, CudaTp<1> >;
      extern template class AmIteratorGrid<2, CudaTp<2> >;
      extern template class AmIteratorGrid<3, CudaTp<3> >;
   }

} 
#endif
