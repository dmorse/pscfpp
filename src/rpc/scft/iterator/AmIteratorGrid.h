#ifndef RPC_AM_ITERATOR_GRID_H
#define RPC_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.h>  // direct base class 
#include <rpc/system/Types.h>                 // direct base argument
#include <rpc/scft/iterator/Iterator.h>       // tertiary base class
//#include <pscf/iterator/AmIteratorTmpl.h>     // secondary base class
#include <prdc/field/cpu/FftwDRArray.h>       // indirect base argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class AmIteratorTmpl<Rp::Iterator<1, Rpc::Types<1> >, Prdc::Cpu::FftwDRArray<double> >;
   extern template 
   class AmIteratorTmpl<Rp::Iterator<2, Rpc::Types<2> >, Prdc::Cpu::FftwDRArray<double> >;
   extern template 
   class AmIteratorTmpl<Rp::Iterator<3, Rpc::Types<3> >, Prdc::Cpu::FftwDRArray<double> >;
   namespace Rp {
      extern template class AmIteratorGrid<1, Rpc::Types<1> >;
      extern template class AmIteratorGrid<2, Rpc::Types<2> >;
      extern template class AmIteratorGrid<3, Rpc::Types<3> >;
   } 
}
#endif
