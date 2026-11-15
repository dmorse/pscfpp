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
#include <rp/solvers/Block.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <pscf/chem/Ensemble.h>
#include <pscf/chem/PolymerModel.h>

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

   namespace {

      /*
      * Correlation factors of one thread block with a discrete contour.
      *
      * The MDE solver takes ns - 1 contour steps of length
      * ds = length/(ns - 1), each by Richardson extrapolation of one
      * full and two half pseudo-spectral steps, and integrates the
      * concentration over the ns contour points by Simpson's rule. To
      * first order in a change of the field of a homogeneous state, a
      * step then acts as if the field were applied at its two ends and
      * its midpoint with weights ds/6, 4ds/6 and ds/6. With
      * x = kSq*kuhn*kuhn/6 and r = exp(-x*ds), this gives
      *
      *   eC = sum_j c_j r^j  (Simpson weights c_j of the contour points)
      *   eW = s*(1 - r^(ns-1))/(1 - r),  s = ds*(1 + 4*r^(1/2) + r)/6
      *   d  = 2*s*(length - eC)/(1 - r)
      *
      * Here d is the intra-block correlation (length^2 times the Debye
      * function for a continuous chain), and eC and eW are the end
      * factors of the concentration points and of the field points
      * (both (1 - exp(-x*length))/x for a continuous chain). All three
      * tend to their continuous values as ds -> 0, but d tends to
      * length*ds/3 rather than 2*length/x for x*ds >> 1.
      */
      void threadFactors(double kSq, double length, double kuhn, int ns,
                         double& d, double& eC, double& eW)
      {
         const double x = kSq*kuhn*kuhn/6.0;
         if (x == 0.0) {
            d = length*length;
            eC = length;
            eW = length;
            return;
         }
         const int n = ns - 1;                 // even number of steps
         const double ds = length/double(n);
         const double a = x*ds;
         const double xl = x*length;
         const double r = std::exp(-a);
         const double omr = -std::expm1(-a);        // 1 - r
         const double omr2 = -std::expm1(-2.0*a);   // 1 - r^2
         const double s = ds*(1.0 + 4.0*std::exp(-0.5*a) + r)/6.0;

         // Simpson sum, with geometric sums over odd j and over even
         // j (0 < j < n), written with expm1 to avoid cancellation
         eC = 1.0 + std::exp(-xl);
         eC += 4.0*r*(-std::expm1(-xl))/omr2;
         eC += 2.0*r*r*(-std::expm1(-(xl - 2.0*a)))/omr2;
         eC *= ds/3.0;

         eW = s*(-std::expm1(-xl))/omr;

         // length - eC, summed term by term where it would cancel
         double diff;
         if (xl < 0.1) {
            diff = -std::expm1(-xl);
            for (int j = 1; j < n; ++j) {
               diff += (j % 2 ? 4.0 : 2.0)*(-std::expm1(-a*double(j)));
            }
            diff *= ds/3.0;
         } else {
            diff = length - eC;
         }
         d = 2.0*s*diff/omr;
      }

   }

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
      if (PolymerModel::isThread()) {
         // Response of the MDE solution with a discrete contour
         computeThreadOmegaTotal(correlations);
      } else {
         // The bead model is discrete as it stands
         correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);
      }

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
   * Compute the total correlation for the thread model (on host).
   */
   template <int D, class T>
   void
   IntraCorrelation<D,T>::computeThreadOmegaTotal(Array<RealT>& correlations)
   {
      Mixture<D,T> const & mixture = system().mixture();
      const double vMonomer = mixture.vMonomer();
      UTIL_CHECK(correlations.capacity() == kSize_);
      UTIL_CHECK(Gsq_.capacity() == kSize_);

      for (int j = 0; j < kSize_; ++j) {
         correlations[j] = 0.0;
      }

      // Polymer species: the continuous-chain sum over pairs of blocks,
      // with the factors of threadFactors for the Debye functions and
      // end factors of the blocks.
      DArray<double> d, eC, eW;
      for (int ip = 0; ip < mixture.nPolymer(); ++ip) {
         Polymer<D,T> const & polymer = mixture.polymer(ip);
         Correlation::Polymer<RealT> const & cPolymer
                                     = correlationMixturePtr_->polymer(ip);
         const int nBlock = polymer.nBlock();
         const double prefactor
                    = cPolymer.phi()/(cPolymer.totalLength()*vMonomer);
         d.allocate(nBlock);
         eC.allocate(nBlock);
         eW.allocate(nBlock);
         for (int j = 0; j < kSize_; ++j) {
            const double kSq = Gsq_[j];
            for (int ia = 0; ia < nBlock; ++ia) {
               threadFactors(kSq, polymer.block(ia).length(),
                             polymer.block(ia).kuhn(),
                             polymer.block(ia).ns(),
                             d[ia], eC[ia], eW[ia]);
            }
            double omega = 0.0;
            for (int ia = 0; ia < nBlock; ++ia) {
               omega += d[ia];
               for (int ib = 0; ib < ia; ++ib) {
                  omega += std::exp(-cPolymer.rSq(ia, ib)*kSq/6.0)
                           * (eC[ia]*eW[ib] + eW[ia]*eC[ib]);
               }
            }
            correlations[j] += prefactor*omega;
         }
         d.deallocate();
         eC.deallocate();
         eW.deallocate();
      }

      // Solvent species (point particles)
      for (int i = 0; i < mixture.nSolvent(); ++i) {
         const double c = mixture.solvent(i).phi()*mixture.solvent(i).size()
                          / vMonomer;
         for (int j = 0; j < kSize_; ++j) {
            correlations[j] += c;
         }
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
