/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/cpu/FftwDArray.h>
#include <pscf/backend/CPT.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>
#include <rp/fts/compressor/AmCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rp::Compressor<1,CPT>, FftwDRArray<double> >;
   template class AmIteratorTmpl< Rp::Compressor<2,CPT>, FftwDRArray<double> >; 
   template class AmIteratorTmpl< Rp::Compressor<3,CPT>, FftwDRArray<double> >;
   namespace Rp {
      template class AmCompressor<1,CPT>;
      template class AmCompressor<2,CPT>;
      template class AmCompressor<3,CPT>;
   }
}
