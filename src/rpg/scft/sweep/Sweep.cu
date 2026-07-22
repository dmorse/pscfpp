/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"

#include <rpg/system/System.h>
#include <rpg/scft/iterator/Iterator.h>
#include <rpg/scft/ScftThermo.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <rpg/field/CFields.h>

#include <rp/scft/sweep/Sweep.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template class SweepTmpl< Rp::BasisFieldState<1,CUT> >;
   template class SweepTmpl< Rp::BasisFieldState<2,CUT> >;
   template class SweepTmpl< Rp::BasisFieldState<3,CUT> >;
   namespace Rp {
      template class Sweep<1,CUT>;
      template class Sweep<2,CUT>;
      template class Sweep<3,CUT>;
   }
}
