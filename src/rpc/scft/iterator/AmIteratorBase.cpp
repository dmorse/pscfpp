/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBase.h"
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rpc::Iterator<1>, DRArray<double> >;
   template class AmIteratorTmpl< Rpc::Iterator<2>, DRArray<double> >; 
   template class AmIteratorTmpl< Rpc::Iterator<3>, DRArray<double> >;
}
