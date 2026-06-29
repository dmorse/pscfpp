#ifndef CPC_MIXTURE_TPP
#define CPC_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include "Polymer.h"
#include "Solvent.h"
#include "Block.h"
#include "Propagator.h"
#include <cpc/field/FieldIo.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/complex.h>

#include <cp/solvers/Mixture.tpp>

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;

   template <int D>
   void Mixture<D>::eqS(FieldT& A, double c) const
   {
      const int nx = mesh().size();
      UTIL_CHECK(nx == A.capacity());
      for (int i = 0; i < nx; ++i) {
         assign(A[i], c);
      }
   }

   template <int D>
   void Mixture<D>::addEqV(FieldT& A, FieldT const & B) const
   {
      const int nx = mesh().size();
      UTIL_CHECK(nx == A.capacity());
      UTIL_CHECK(nx == B.capacity());
      for (int i = 0; i < nx; ++i) {
         assign(A[i], B[i]);
      }
   }

   /*
   * Allocate memory for all blocks.
   */
   template <int D>
   void Mixture<D>::allocateBlocks()
   {
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).allocate(ds());
         }
      }
   }

} // namespace Cpc
} // namespace Pscf
#endif
