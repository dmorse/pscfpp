#ifndef RPC_MIXTURE_TPP
#define RPC_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <rpc/field/FieldIo.h>
#include <prdc/solvers/MixtureReal.tpp>
#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>

namespace Pscf {
namespace Rpc {

   using namespace Prdc;

   template <int D>
   void Mixture<D>::eqS(FieldT& A, double c) const
   {
      const int nx = mesh().size();
      UTIL_CHECK(nx == A.capacity());
      for (int i = 0; i < nx; ++i) {
         A[i] = c;
      }
   }

   template <int D>
   void Mixture<D>::addEqV(FieldT& A, FieldT const & B) const
   {
      const int nx = mesh().size();
      UTIL_CHECK(nx == A.capacity());
      UTIL_CHECK(nx == B.capacity());
      for (int i = 0; i < nx; ++i) {
         A[i] += B[i];
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

} // namespace Rpc
} // namespace Pscf
#endif
