/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/System.h>
#include <rp/scft/iterator/Iterator.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/backends/CPT.h>
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
