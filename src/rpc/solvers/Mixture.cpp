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
   void Mixture<D, Cpp<D> >::allocateBlocks()
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
   class MixtureTmpl< Rp::Polymer<1, Cpp<1> >, Rp::Solvent<1, Cpp<1> > >;
   template 
   class MixtureTmpl< Rp::Polymer<2, Cpp<2> >, Rp::Solvent<2, Cpp<2> > >;
   template 
   class MixtureTmpl< Rp::Polymer<3, Cpp<3> >, Rp::Solvent<3, Cpp<3> > >;
   namespace Rp { 
      template class MixtureBase<1, Cpp<1> >;
      template class MixtureBase<2, Cpp<2> >;
      template class MixtureBase<3, Cpp<3> >;
      template class Mixture<1, Cpp<1> >;
      template class Mixture<2, Cpp<2> >;
      template class Mixture<3, Cpp<3> >;
   }
}
