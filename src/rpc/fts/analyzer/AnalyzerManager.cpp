/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/analyzer/AnalyzerManager.h>
#include <rpc/fts/analyzer/AnalyzerFactory.h>

#include <rp/fts/analyzer/AnalyzerManager.tpp> // class implementation

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class AnalyzerManager<1, Rpc::Types<1> >;
      template class AnalyzerManager<2, Rpc::Types<2> >;
      template class AnalyzerManager<3, Rpc::Types<3> >;
   }
}
