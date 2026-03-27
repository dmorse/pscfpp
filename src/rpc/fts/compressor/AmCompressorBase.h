#ifndef RPC_AM_COMP_BASE_H
#define RPC_AM_COMP_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* Instantiatied class template specializations used as base classes 
* for specializations of the class templates Rpc::AmCompressor and 
* Rpc::LrAmCompressor.
*/

#include <pscf/iterator/AmIteratorTmpl.h>   // base class template
#include <rpc/fts/compressor/Compressor.h>  // base class argument
#include <util/containers/DArray.h>         // base class argument

namespace Pscf {
   extern template 
   class AmIteratorTmpl< Rpc::Compressor<1>, DArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpc::Compressor<2>, DArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpc::Compressor<3>, DArray<double> >;
}
#endif
