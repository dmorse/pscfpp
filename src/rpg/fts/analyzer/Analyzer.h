#ifndef RPG_ANALYZER_H
#define RPG_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>   // class template
#include <pscf/cuda/Cuda.h>           // class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1, CudaTp<1> >;
      extern template class Analyzer<2, CudaTp<2> >;
      extern template class Analyzer<3, CudaTp<3> >;
   }
}
#endif
