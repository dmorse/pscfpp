/*
* PSCF - MixtureModifier Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureModifier.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/solvers/Block.h>

#include <rp/solvers/MixtureModifier.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MixtureModifier<1,CUT>;
      template class MixtureModifier<2,CUT>;
      template class MixtureModifier<3,CUT>;
   }
}
