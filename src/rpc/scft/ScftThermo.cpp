/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"                     // class header
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <rpc/field/Mask.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/ScftThermo.tpp>           // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1, CppTp<1> >;
      template class ScftThermo<2, CppTp<2> >;
      template class ScftThermo<3, CppTp<3> >;
   }
}
