#ifndef RP_INTRACORRELATION_TPP
#define RP_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/correlation/Mixture.h>

#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   IntraCorrelation<D,T>::IntraCorrelation(
                                     typename T::System const & system)
    : systemPtr_(&system),
      correlationMixturePtr_(nullptr),
      kSize_(-1)
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
   * Compute k-space array of intramolecular correlation functions.
   */
   template <int D, class T>
   void
   IntraCorrelation<D,T>::computeOmegaTotal(Array<RealT>& correlations)
   {
      computeMeshProperties();

      // Compute total intramolecular correlation function
      UTIL_CHECK(correlations.capacity() == kSize_);
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);

   }

   /*
   * Compute k-grid mesh dimensions and Gsq.
   */
   template <int D, class T>
   void IntraCorrelation<D,T>::computeMeshProperties()
   {
      // Local copies of domain properties
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Compute k-space mesh dimensions kMeshDimensions_ and size Size_
      typename T::FFT:computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Check allocation of Gsq_ (k-space array of square wavenumbers)
      if (!Gsq_.isAllocated()) {
         Gsq_.allocate(kSize_);
      }
      UTIL_CHECK(Gsq_.capacity() == kSize_);

      // Compute Gsq_
      IntVec<D> G, Gmin;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, dimensions, unitCell);
         Gsq_[iter.rank()] = unitCell.ksq(Gmin);
      }
   }

}
}
#endif
