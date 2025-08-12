#ifndef PRDC_UNIT_CELL_BASE_TPP
#define PRDC_UNIT_CELL_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCellBase.h"

#include <util/signal/Signal.h>
#include <util/math/Constants.h>
#include <util/containers/FSArray.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   UnitCellBase<D>::UnitCellBase()
    : nParameter_(0),
      isInitialized_(false),
      signalPtr_(nullptr)
   {  initializeToZero(); }

   /*
   * Destructor.
   */
   template <int D>
   UnitCellBase<D>::~UnitCellBase()
   {}

   /*
   * Set all the unit cell parameters.
   */
   template <int D>
   void 
   UnitCellBase<D>::setParameters(FSArray<double, 6> const& parameters)
   {
      UTIL_CHECK(parameters.size() == nParameter_);
      isInitialized_ = false;
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i] = parameters[i];
      }
      setLattice();
   }

   /*
   * Get all the unit cell parameters.
   */
   template <int D>
   FSArray<double, 6> UnitCellBase<D>::parameters() const
   {
      FSArray<double, 6> parameters;
      for (int i = 0; i < nParameter_; ++i) {
         parameters.append(parameters_[i]);
      }
      return parameters;
   }
   
   /*
   * Get square magnitude of reciprocal basis vector.
   */
   template <int D>
   double UnitCellBase<D>::ksq(IntVec<D> const & k) const
   {
      RealVec<D> g(0.0);
      RealVec<D> p;
      for (int i = 0; i < D; ++i) {
         p.multiply(kBasis_[i], k[i]);
         g += p;
      }
      double value = 0.0;
      for (int i = 0; i < D; ++i) {
         value += g[i]*g[i];
      }
      return value;
   }

   /*
   * Get magnitude of derivative of square of reciprocal basis vector.
   */
   template <int D>
   double UnitCellBase<D>::dksq(IntVec<D> const & vec, int n) const
   {
      double element = 0.0;
      double value = 0.0;

      for (int p = 0; p < D; ++p){
         for (int q = 0; q < D; ++q){
            element = dkkBasis(n, p, q);
            value += vec[p]*vec[q]*element;
         }
      }

      return value;
   }

   // Signal

   template <int D>
   void UnitCellBase<D>::setSignal(Signal<void>& signal)
   {  signalPtr_ = &signal; }

   template <int D>
   bool UnitCellBase<D>::hasSignal() const
   {  return bool(signalPtr_); }

   template <int D>
   Signal<void>& UnitCellBase<D>::signal()
   {  
      UTIL_CHECK(signalPtr_);
      return *signalPtr_;
   }

   // Protected member functions

   /*
   * Initialize internal arrays to zero.
   */
   template <int D>
   void UnitCellBase<D>::initializeToZero()
   {
      // Initialize all elements to zero
      int i, j, k;
      for (i = 0; i < D; ++i) {
         for (j = 0; j < D; ++j) {
            rBasis_[i][j] = 0.0;
            kBasis_[i][j] = 0.0;
         }
      }
      for (k = 0; k < 6; ++k){
         for (i = 0; i < D; ++i) {
            for (j = 0; j < D; ++j) {
               drBasis_[k](i,j) = 0.0;
               dkBasis_[k](i,j) = 0.0;
               drrBasis_[k](i,j) = 0.0;
               dkkBasis_[k](i,j) = 0.0;
            }
         }
      }
   }

   /*
   * Compute quantities involving derivatives.
   */
   template <int D>
   void UnitCellBase<D>::computeDerivatives()
   {
      // Compute dkBasis
      int p, q, r, s, t;
      for (p = 0; p < nParameter_; ++p) {
         for (q = 0; q < D; ++q) {
            for (r = 0; r < D; ++r) {

               // Loop over free indices s, t
               for (s = 0; s < D; ++s) {
                  for (t = 0; t < D; ++t) {
                     dkBasis_[p](q,r)
                       -= kBasis_[q][s]*drBasis_[p](t,s)*kBasis_[t][r];
                  }
               }
               dkBasis_[p](q,r) /= 2.0*Constants::Pi;

            }
         }
      }

      // Compute drrBasis and dkkBasis
      for (p = 0; p < nParameter_; ++p) {
         for (q = 0; q < D; ++q) {
            for (r = 0; r < D; ++r) {
               for (s = 0; s < D; ++s) {
                  drrBasis_[p](q,r) += rBasis_[q][s]*drBasis_[p](r,s);
                  drrBasis_[p](q,r) += rBasis_[r][s]*drBasis_[p](q,s);
                  dkkBasis_[p](q,r) += kBasis_[q][s]*dkBasis_[p](r,s);
                  dkkBasis_[p](q,r) += kBasis_[r][s]*dkBasis_[p](q,s);

               }
            }
         }
      }

   }

   /*
   * Set all lattice parameters.
   */
   template <int D>
   void UnitCellBase<D>::setLattice()
   {
      initializeToZero();
      setBasis();
      computeDerivatives();
      isInitialized_ = true;
      if (hasSignal()) {
         signalPtr_->notify();
      }
   }

}
}
#endif
