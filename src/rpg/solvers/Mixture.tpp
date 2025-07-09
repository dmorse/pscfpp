#ifndef RPG_MIXTURE_TPP
#define RPG_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <prdc/solvers/MixtureReal.tpp>
#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>


namespace Pscf {
namespace Rpg {

   using namespace Prdc;
   using namespace Prdc::Cuda;

   /*
   * Constructor
   */
   template <int D>
   Mixture<D>::Mixture()
    : MixtureRealT(),
      useBatchedFFT_(true)
   {}

   /*
   * Read all parameters and initialize.
   */
   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureRealT::readParameters(in);

      // Optionally read useBatchedFFT boolean
      useBatchedFFT_ = true;
      ParamComposite::readOptional(in, "useBatchedFFT", useBatchedFFT_);
   }

   /*
   * Set all elements of a field to a single scalar: A[i] = c.
   */
   template <int D>
   void Mixture<D>::eqS(FieldT& A, double c) const
   {
      const int nx = mesh().size();
      UTIL_CHECK(nx == A.capacity());
      VecOp::eqS(A,c);
   }

   /*
   * Compound addition-assignment of two fields: A[i] += B[i]
   */
   template <int D>
   void Mixture<D>::addEqV(FieldT& A, FieldT const & B) const
   {
      const int nx = mesh().size();
      UTIL_CHECK(nx == A.capacity());
      UTIL_CHECK(nx == B.capacity());
      VecOp::addEqV(A, B);
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
            polymer(i).block(j).allocate(ds(), useBatchedFFT_);
         }
      }
   }

} // namespace Rpg
} // namespace Pscf
#endif
