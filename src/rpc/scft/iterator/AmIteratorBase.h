#ifndef RPC_AM_ITERATOR_BASE_H
#define RPC_AM_ITERATOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* Declarations of explicit instantations of AmIteratorTmpl used as base 
* classes for class templates Rpc::AmIteratorBasis and Rpc::AmIteratorGrid.
*/

#include <pscf/iterator/AmIteratorTmpl.h>   // base class template
#include <rpc/scft/iterator/Iterator.h>     // base class argument
#include <util/containers/DRArray.h>        // base class argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class AmIteratorTmpl< Rpc::Iterator<1>, DRArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpc::Iterator<2>, DRArray<double> >;
   extern template 
   class AmIteratorTmpl< Rpc::Iterator<3>, DRArray<double> >;
}
#endif
