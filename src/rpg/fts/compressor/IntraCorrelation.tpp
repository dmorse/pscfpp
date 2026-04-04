#ifndef RPG_INTRACORRELATION_TPP
#define RPG_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/correlation/Mixture.h>

#include <util/global.h>

namespace Pscf {
namespace Rpg{

   using namespace Util;
   using namespace Prdc; 
   using namespace Prdc::Cuda; 

   /*
   * Constructor.
   */
   template <int D>
   IntraCorrelation<D>::IntraCorrelation(System<D> const & system)
    : systemPtr_(&system),
      correlationMixturePtr_(nullptr),
      kSize_(-1)
   {
      correlationMixturePtr_ 
          = new Correlation::Mixture<cudaReal>(system.mixture());
   }

   /*
   * Destructor.
   */
   template <int D>
   IntraCorrelation<D>::~IntraCorrelation()
   {  delete correlationMixturePtr_; }

   template<int D>
   void
   IntraCorrelation<D>::computeOmegaTotal(RField<D>& correlations)
   {
      computeMeshProperties();

      // Check allocation of HostDArray<double> correlations_
      if (!correlations_.isAllocated()) {
         correlations_.allocate(kSize_);
      }
      UTIL_CHECK(correlations_.capacity() == kSize_);

      // Compute array of correlation function values on CPU
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations_);

      // Copy host array "correlations_" to device array "correlations"
      UTIL_CHECK(correlations.capacity() == kSize_);
      correlations = correlations_;
   }

   /*
   * Compute k-grid mesh dimensions and Gsq.
   */
   template<int D>
   void IntraCorrelation<D>::computeMeshProperties()
   {
      // Local copies of domain properties
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Compute k-space mesh dimensions kMeshDimensions_ and size kSize_
      FFT<D>::computeKMesh(dimensions, kMeshDimensions_, kSize_);

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
