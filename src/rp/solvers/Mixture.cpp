/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOp.h>
#include <rp/solvers/Mixture.tpp>

// Explicit instantiation definitions
namespace Pscf {
   template 
   class MixtureTmpl< Rp::Polymer<1,CPT>, Rp::Solvent<1,CPT> >;
   template 
   class MixtureTmpl< Rp::Polymer<2,CPT>, Rp::Solvent<2,CPT> >;
   template 
   class MixtureTmpl< Rp::Polymer<3,CPT>, Rp::Solvent<3,CPT> >;
   namespace Rp { 
      template class Mixture<1,CPT>;
      template class Mixture<2,CPT>;
      template class Mixture<3,CPT>;
   }
}
