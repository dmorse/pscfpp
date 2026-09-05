/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <rp/fts/analyzer/ConcentrationWriter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationWriter<1,CPT>;
      template class ConcentrationWriter<2,CPT>;
      template class ConcentrationWriter<3,CPT>;
   }
}
