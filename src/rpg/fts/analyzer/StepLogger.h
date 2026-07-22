#ifndef RPG_STEP_LOGGER_H
#define RPG_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/StepLogger.h>  // base class template
#include <pscf/backends/CUT.h>            // base class argument
#include <rpg/fts/analyzer/Analyzer.h>   // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class StepLogger<1,CUT>;
      extern template class StepLogger<2,CUT>;
      extern template class StepLogger<3,CUT>;
   }
}
#endif
