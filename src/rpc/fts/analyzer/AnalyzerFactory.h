#ifndef RPC_ANALYZER_FACTORY_H
#define RPC_ANALYZER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AnalyzerFactory.h>
#include <rpc/system/Types.h>

// Explicit instantiation declarations
namespace Pscf {
namespace Rp {
   extern template class AnalyzerFactory<1, Rpc::Types<1> >;
   extern template class AnalyzerFactory<2, Rpc::Types<2> >;
   extern template class AnalyzerFactory<3, Rpc::Types<3> >;
}
}
#endif
