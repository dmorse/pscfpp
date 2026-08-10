/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AnalyzerManager.tpp> // class implementation

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class AnalyzerManager<1,CUT>;
      template class AnalyzerManager<2,CUT>;
      template class AnalyzerManager<3,CUT>;
   }
}
