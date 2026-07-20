#ifndef RPG_ANALYZER_FACTORY_H
#define RPG_ANALYZER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AnalyzerFactory.h>
#include <pscf/cuda/CudaTp.h>

// Explicit instantiation declarations
namespace Pscf {
namespace Rp {
   extern template class AnalyzerFactory<1, CudaTp<1> >;
   extern template class AnalyzerFactory<2, CudaTp<2> >;
   extern template class AnalyzerFactory<3, CudaTp<3> >;
}
}
#endif
