/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
#include "HamiltonianAnalyzer.h"                    // header
#include <rp/fts/simulator/Simulator.h>
#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#endif

#include <rp/fts/analyzer/HamiltonianAnalyzer.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class HamiltonianAnalyzer<1,CPT>;
      template class HamiltonianAnalyzer<2,CPT>;
      template class HamiltonianAnalyzer<3,CPT>;
   }
}
