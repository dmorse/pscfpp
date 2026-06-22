/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"               // class header

#include "Polymer.h"
#include "Solvent.h"
#include "Block.h"
#include "Propagator.h"
#include <rpg/field/FieldIo.h>
#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/VecOp.h>

#include <rp/solvers/Mixture.tpp>  // base class template implementation

namespace Pscf {
namespace Rpg {

   using namespace Prdc;

   /*
   * Constructor
   */
   template <int D>
   Mixture<D>::Mixture()
    : Rp::Mixture<D, Types<D> >(),
      useBatchedFFT_(true)
   {}

   /*
   * Read all parameters and initialize.
   */
   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      RpMixtureT::readParameters(in);

      // Optionally read useBatchedFFT boolean
      useBatchedFFT_ = true;
      ParamComposite::readOptional(in, "useBatchedFFT", useBatchedFFT_);
   }

   /*
   * Allocate memory for all blocks.
   */
   template <int D>
   void Mixture<D>::allocateBlocks()
   {
      const double ds = RpMixtureT::ds();
      const int np = Composition<cudaReal>::nPolymer();
      int i, j;
      for (i = 0; i < np; ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).allocate(ds, useBatchedFFT_);
         }
      }
   }

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   template 
   class MixtureTmpl< Rp::Polymer<1, Rpg::Types<1> >, Rp::Solvent<1, Rpg::Types<1> > >;
   template 
   class MixtureTmpl< Rp::Polymer<2, Rpg::Types<2> >, Rp::Solvent<2, Rpg::Types<2> > >;
   template 
   class MixtureTmpl< Rp::Polymer<3, Rpg::Types<3> >, Rp::Solvent<3, Rpg::Types<3> > >;
   namespace Rp {
      template class Mixture<1, Rpg::Types<1> >;
      template class Mixture<2, Rpg::Types<2> >;
      template class Mixture<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}
