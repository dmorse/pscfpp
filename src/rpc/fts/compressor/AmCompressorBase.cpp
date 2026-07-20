/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressorBase.h"

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rp::Compressor<1, CppTp<1> >, FftwDRArray<double> >;
   template class AmIteratorTmpl< Rp::Compressor<2, CppTp<2> >, FftwDRArray<double> >; 
   template class AmIteratorTmpl< Rp::Compressor<3, CppTp<3> >, FftwDRArray<double> >;
}
