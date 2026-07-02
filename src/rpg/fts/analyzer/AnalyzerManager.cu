/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerManager.h"
#include "AnalyzerFactory.h"
#include <rp/fts/analyzer/AnalyzerManager.tpp> // class implementation

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class AnalyzerManager<1, Rpg::Types<1> >;
      template class AnalyzerManager<2, Rpg::Types<2> >;
      template class AnalyzerManager<3, Rpg::Types<3> >;
   }
}
