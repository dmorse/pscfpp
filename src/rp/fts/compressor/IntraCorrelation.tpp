#ifndef RP_INTRACORRELATION_TPP
#define RP_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <pscf/chem/Ensemble.h>

#include <pscf/field/FFT.h>
#include <pscf/field/RField.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/correlation/Mixture.h>
#include <pscf/correlation/Polymer.h>

#include <util/global.h>

#include <cmath>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   IntraCorrelation<D,T>::IntraCorrelation(
                                     System<D,T> const & system)
    : systemPtr_(&system),
      correlationMixturePtr_(nullptr),
      kSize_(-1),
      isCurrent_(false)
   {
      correlationMixturePtr_
          = new Pscf::Correlation::Mixture<RealT>(system.mixture());
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   IntraCorrelation<D,T>::~IntraCorrelation()
   {  delete correlationMixturePtr_; }

   /*
   * Compute k-space array of intramolecular correlation functions (on device).
   */
   template<int D, class T>
   void
   IntraCorrelation<D,T>::computeOmegaTotal(RField<D,T>& correlations)
   {
      getMeshDimensions();
      int nk = kSize();
      UTIL_CHECK(correlations.isAllocated());
      UTIL_CHECK(correlations.capacity() == nk);

      correlations_h_.associate(correlations);
      computeOmegaTotalArray(correlations_h_);
      correlations = correlations_h_;

      // Release the association, so that the correlations can be
      // computed again
      correlations_h_.dissociate();
   }

   /*
   * Compute the correlation function if it may have changed.
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::update(RField<D,T>& correlations)
   {
      // Volume fractions of open species, which depend on the fields
      Mixture<D,T> const & mixture = system().mixture();
      std::vector<double> openPhi;
      for (int i = 0; i < mixture.nPolymer(); ++i) {
         if (mixture.polymer(i).ensemble() == Ensemble::Open) {
            openPhi.push_back(mixture.polymer(i).phi());
         }
      }
      for (int i = 0; i < mixture.nSolvent(); ++i) {
         if (mixture.solvent(i).ensemble() == Ensemble::Open) {
            openPhi.push_back(mixture.solvent(i).phi());
         }
      }

      // Recompute if one of them changed by more than 0.1 percent, or
      // if a stored value is not finite (from a diverging iteration)
      if (isCurrent_) {
         UTIL_CHECK(openPhi.size() == openPhi_.size());
         for (size_t i = 0; i < openPhi.size(); ++i) {
            double change = std::abs(openPhi[i] - openPhi_[i]);
            if (!std::isfinite(openPhi_[i]) || change > 1.0E-3*openPhi_[i]) {
               isCurrent_ = false;
            }
         }
      }

      if (!isCurrent_) {
         computeOmegaTotal(correlations);
         openPhi_ = openPhi;
         isCurrent_ = true;
      }
   }

   /*
   * Mark the correlation function for recomputation.
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::clearCache()
   {  isCurrent_ = false; }

   /*
   * Compute k-space array of intramolecular correlation functions (on host).
   */
   template <int D, class T>
   void
   IntraCorrelation<D,T>::computeOmegaTotalArray(Array<RealT>& correlations)
   {
      getMeshDimensions();
      UTIL_CHECK(correlations.capacity() == kSize_);
      computeGsq();
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);

      // The LR compressors use the correlation as the response of the
      // total concentration to the pressure-like field. A uniform change
      // of this field does not change the concentration of a closed
      // species, whose amount is fixed. In a mixed ensemble, the k = 0
      // value (rank 0) is thus the sum of the values phi*N/v of the
      // open species alone. (In a closed system the k = 0 residual is
      // zero, and in an open one every species responds.)
      Mixture<D,T> const & mixture = system().mixture();
      const double vMonomer = mixture.vMonomer();
      bool hasOpen = false;
      bool hasClosed = false;
      double omegaOpen = 0.0;
      for (int i = 0; i < mixture.nPolymer(); ++i) {
         if (mixture.polymer(i).ensemble() == Ensemble::Open) {
            Correlation::Polymer<RealT> const & polymer
                                     = correlationMixturePtr_->polymer(i);
            omegaOpen += polymer.phi()*polymer.totalLength()/vMonomer;
            hasOpen = true;
         } else {
            hasClosed = true;
         }
      }
      for (int i = 0; i < mixture.nSolvent(); ++i) {
         if (mixture.solvent(i).ensemble() == Ensemble::Open) {
            omegaOpen += mixture.solvent(i).phi()*mixture.solvent(i).size()
                         / vMonomer;
            hasOpen = true;
         } else {
            hasClosed = true;
         }
      }
      if (hasOpen && hasClosed) {
         correlations[0] = omegaOpen;
      }
   }

   /*
   * Compute k-grid mesh dimensions and Gsq.
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::getMeshDimensions()
   {
      meshDimensions_ = system().domain().mesh().dimensions();

      // Compute k-space mesh dimensions kMeshDimensions_ and size Size_
      FFT<D,T>::computeKMesh(meshDimensions_, kMeshDimensions_, kSize_);
   }

   /*
   * Construct array of squared wavevector values (on host).
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::computeGsq()
   {
      // Get mesh dimensions if not done previously
      if (kSize_ <= 0) {
         getMeshDimensions();
      }

      // Check allocation of Gsq_ (k-space array of square wavenumbers)
      if (!Gsq_.isAllocated()) {
         Gsq_.allocate(kSize_);
      }
      UTIL_CHECK(Gsq_.capacity() == kSize_);

      // Keep Gsq_ if the unit cell parameters have not changed
      UnitCell<D> const & unitCell = system().domain().unitCell();
      FSArray<double, 6> const parameters = unitCell.parameters();
      bool isCurrent = (parameters.size() == gsqParameters_.size());
      for (int i = 0; isCurrent && i < parameters.size(); ++i) {
         if (parameters[i] != gsqParameters_[i]) isCurrent = false;
      }
      if (isCurrent) return;
      gsqParameters_ = parameters;

      // Compute Gsq_
      IntVec<D> G, Gmin;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, meshDimensions_, unitCell);
         Gsq_[iter.rank()] = unitCell.ksq(Gmin);
      }
   }

}
}
#endif
