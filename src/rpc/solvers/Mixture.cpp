/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/solvers/Block.h>
#include <rpc/solvers/Propagator.h>
#include <rpc/field/FieldIo.h>

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>

#include <pscf/cpu/VecOp.h>

#include <rp/solvers/MixtureBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Prdc;

   /*
   * Allocate memory for all blocks.
   */
   template <int D>
   void Mixture<D,CPT>::allocateBlocks()
   {
      const double ds = RpMixtureT::ds();
      const int np = CompositionT::nPolymer();
      int i, j;
      for (i = 0; i < np; ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).allocate(ds);
         }
      }
   }

} 
} 

// Explicit instantiation definitions
namespace Pscf {
   template 
   class MixtureTmpl< Rp::Polymer<1,CPT>, Rp::Solvent<1,CPT> >;
   template 
   class MixtureTmpl< Rp::Polymer<2,CPT>, Rp::Solvent<2,CPT> >;
   template 
   class MixtureTmpl< Rp::Polymer<3,CPT>, Rp::Solvent<3,CPT> >;
   namespace Rp { 
      template class MixtureBase<1,CPT>;
      template class MixtureBase<2,CPT>;
      template class MixtureBase<3,CPT>;
      template class Mixture<1,CPT>;
      template class Mixture<2,CPT>;
      template class Mixture<3,CPT>;
   }
}
