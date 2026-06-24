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
   template class SweepTmpl< Rp::BasisFieldState<1, Rpg::Types<1> > >;
   template class SweepTmpl< Rp::BasisFieldState<2, Rpg::Types<2> > >;
   template class SweepTmpl< Rp::BasisFieldState<3, Rpg::Types<3> > >;
   namespace Rp {
      template class Sweep<1, Rpg::Types<1> >;
      template class Sweep<2, Rpg::Types<2> >;
      template class Sweep<3, Rpg::Types<3> >;
   }
}
