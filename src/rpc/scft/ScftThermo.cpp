/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/ScftThermo.tpp>           // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1,CPT>;
      template class ScftThermo<2,CPT>;
      template class ScftThermo<3,CPT>;
   }
}
