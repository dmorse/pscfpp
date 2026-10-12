/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecOp.h>
#include <rp/solvers/Mixture.tpp> 

// Explicit instantiation definitions
namespace Pscf {
   template 
   class MixtureTmpl< Rp::Polymer<1,CUT>, Rp::Solvent<1,CUT> >;
   template 
   class MixtureTmpl< Rp::Polymer<2,CUT>, Rp::Solvent<2,CUT> >;
   template 
   class MixtureTmpl< Rp::Polymer<3,CUT>, Rp::Solvent<3,CUT> >;
   namespace Rp {
      template class Mixture<1,CUT>;
      template class Mixture<2,CUT>;
      template class Mixture<3,CUT>;
   }
}
