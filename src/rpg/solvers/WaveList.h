#ifndef RPG_WAVE_LIST_H
#define RPG_WAVE_LIST_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RFieldDft.h>
#include <prdc/cuda/RField.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>
#include <util/containers/DMatrix.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Container for wavevector data.
   */
   template <int D>
   class WaveList
   {
   public:

      /**
      * Constructor.
      */
      WaveList();

      /**
      * Destructor
      */
      ~WaveList();

      /**
      * Allocate memory for all arrays. 
      *
      * \param mesh  spatial discretization mesh (input)
      * \param unitCell  crystallographic unit cell (input)
      */
      void allocate(Mesh<D> const & mesh, 
                    UnitCell<D> const & unitCell);

      /**
      * Compute minimum images of wavevectors.
      *
      * This is only done once, at beginning, using mesh and the initial 
      * unit cell. This can also be called in the readParameters function.
      *
      * \param mesh  spatial discretization Mesh<D> object
      * \param unitCell  crystallographic UnitCell<D> object
      */
      void computeMinimumImages(Mesh<D> const & mesh, 
                                UnitCell<D> const & unitCell);

      /**
      * Compute square norm |k|^2 for all wavevectors.
      *
      * Called once per iteration with unit cell relaxation.
      * Implementation can copy geometric data (rBasis, kBasis, etc.)
      * into local data structures and implement the actual 
      * calculation of kSq or dKSq for each wavevector on the GPU.
      *
      * \param unitCell crystallographic UnitCell<D>
      */
      void computeKSq(UnitCell<D> const & unitCell);

      /**
      * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
      *
      * Called once per iteration with unit cell relaxation.
      *
      * \param unitCell crystallographic UnitCell<D>
      */
      void computedKSq(UnitCell<D> const & unitCell);

      /**
      * Get the minimum image vector for a specified wavevector.
      *
      * \param i index for wavevector
      */
      const IntVec<D>& minImage(int i) const;

      /**
      * Get the kSq array on the host by reference.
      */
      HostDArray<cudaReal> const & kSq() const;

      /**
      * Get the dkSq array on the device by reference.
      */
      DeviceArray<cudaReal> const & dkSq() const;

      /**
      * Get size of k-grid (number of wavewavectors).
      */
      int kSize() const;

      /**
      *  Has memory been allocated for arrays?
      */ 
      bool isAllocated() const
      {  return isAllocated_; }

      /**
      *  Have minimum images been computed?
      */ 
      bool hasMinimumImages() const
      {  return hasMinimumImages_; }

   private:

      // Array containing precomputed minimum images
      DeviceArray<int> minImage_;
      DArray< IntVec<D> > minImage_h_;

      // Array containing values of kSq_, stored on the host
      HostDArray<cudaReal> kSq_h_;

      // Array containing values of dkSq_ stored on the device
      DeviceArray<cudaReal>  dkSq_;

      DeviceArray<cudaReal> dkkBasis_;
      HostDArray<cudaReal> dkkBasis_h_;

      DeviceArray<int> partnerIdTable_;
      HostDArray<int> partnerIdTable_h_;

      DeviceArray<int> selfIdTable_;
      HostDArray<int> selfIdTable_h_;

      DeviceArray<bool> implicit_;
      HostDArray<bool> implicit_h_;

      IntVec<D> dimensions_;
      int kSize_;
      int rSize_;
      int nParams_;

      bool isAllocated_;

      bool hasMinimumImages_;

   };

   template <int D>
   inline const IntVec<D>& WaveList<D>::minImage(int i) const
   {  
      UTIL_CHECK(hasMinimumImages_);
      return minImage_h_[i]; 
   }

   template <int D>
   inline HostDArray<cudaReal> const & WaveList<D>::kSq() const
   {  return kSq_h_; }

   template <int D>
   inline DeviceArray<cudaReal> const & WaveList<D>::dkSq() const
   {  return dkSq_; }

   template <int D>
   inline int WaveList<D>::kSize() const
   { return kSize_; }

   #ifndef RPG_WAVE_LIST_TPP
   // Suppress implicit instantiation
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;
   #endif

}
}
//#include "WaveList.tpp"
#endif
