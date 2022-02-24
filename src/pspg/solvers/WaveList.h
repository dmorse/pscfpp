#ifndef PSPG_WAVE_LIST_H
#define PSPG_WAVE_LIST_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>
#include <pspg/field/RDFieldDft.h>
#include <pspg/field/RDField.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>
#include <util/containers/DMatrix.h>

namespace Pscf { 
namespace Pspg
{ 

   using namespace Util;

   template <int D>
   class WaveList{
   public:
      
      WaveList();
      ~WaveList();

      // Allocate memory for all arrays. Call in readParameters.
      void allocate(Mesh<D>& mesh, UnitCell<D>& unitCell);
         
      // Compute minimum images - only done once, at beginning,
      // using mesh and the initial unit cell. This can also be
      // called in the readParameters function.
      void computeMinimumImages(Mesh<D>& mesh, UnitCell<D>& unitCell);

      // These are called once per iteration with unit cell relaxation
      // Implementation can copy geometric data (rBasis, kBasis, etc.)
      // into local data structures and implement the actual calculation
      // of kSq or dKSq for each wavevector on the GPU. 
      
      void computeKSq(UnitCell<D>& unitCell);
      void computedKSq(UnitCell<D>& unitCell);

      // Accessors: Return pointer to arrays, which are allocated in 
      // the allocate function(). I assume this is the form needed
      // to allow fast use in GPU kernels. 
      
      const IntVec<D>& minImage(int i) const;
      cudaReal* kSq() const;
      cudaReal* dkSq() const;
      int kSize() const;

   private:

      // Bare C array holding precomputed minimum images
      //int* minImage_;
      int* minImage_d;
      // Bare C array holding values of kSq_
      cudaReal*  kSq_;

      // Bare C array holding values of dkSq_
      cudaReal*  dkSq_;

      cudaReal* dkkBasis_d;
      cudaReal* dkkBasis;

      int* partnerIdTable;
      int* partnerIdTable_d;

      int* selfIdTable;
      int* selfIdTable_d;

      bool* implicit;
      bool* implicit_d;

      IntVec<D> dimensions_;
      int kSize_;
      int rSize_;
      int nParams_;

      DArray< IntVec<D> > minImage_;

      bool deviceIsAllocated_;
   };

   template <int D>
   inline const IntVec<D>& WaveList<D>::minImage(int i) const
   { return minImage_[i]; }

   template <int D>
   inline cudaReal* WaveList<D>::kSq() const
   { return kSq_; }

   template <int D>
   inline cudaReal* WaveList<D>::dkSq() const
   { return dkSq_; }

   template <int D>
   inline int WaveList<D>::kSize() const
   { return kSize_; }

   #ifndef PSPG_WAVE_LIST_TPP
   // Suppress implicit instantiation
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;
   #endif

}
}
//#include "WaveList.tpp"
#endif
