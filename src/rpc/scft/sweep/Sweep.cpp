/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"

#include <rpc/system/System.h>
#include <rpc/scft/iterator/Iterator.h>
#include <rpc/scft/ScftThermo.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <rp/scft/sweep/Sweep.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template class SweepTmpl< Rp::BasisFieldState<1,CPT> >;
   template class SweepTmpl< Rp::BasisFieldState<2,CPT> >;
   template class SweepTmpl< Rp::BasisFieldState<3,CPT> >;
   namespace Rp {
      template class Sweep<1,CPT>;
      template class Sweep<2,CPT>;
      template class Sweep<3,CPT>;
   }
}
