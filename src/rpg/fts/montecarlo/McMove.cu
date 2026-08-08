/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"

#include <rpg/fts/montecarlo/McSimulator.h>
#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <rp/fts/montecarlo/McMove.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMove<1,CUT>;
      template class McMove<2,CUT>;
      template class McMove<3,CUT>;
   }
}
