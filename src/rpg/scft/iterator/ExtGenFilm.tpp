#ifndef RPG_EXT_GEN_FILM_TPP
#define RPG_EXT_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExtGenFilm.h"
#include "Iterator.h"
#include <rpg/field/FieldIo.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/paramIdConversions.h>
#include <prdc/cuda/resources.h>
#include <cmath>

namespace Pscf {
namespace Rpg
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   // CUDA kernels: 
   // (defined in anonymous namespace, used only in this file)

   namespace {

      /*
      * Compute the derivative of h fields with respect to film thickness.
      *
      * (note: indices are inverted from their RField representation. So,
      * for a 3D system, x, y, and z in an RField become z, y, and x in the 
      * CUDA thread grid, respectively. In 2D, x, y becomes y, x. This 
      * allows for memory coalescing; see ThreadMesh::setConfig 
      * documentation for more information.)
      * 
      * \param hDerivs  array of fields containing the h derivatives (output)
      * \param mask  field containing the mask
      * \param chiBottom  array of chiBottom values for each species
      * \param chiTop  array of chiTop values for each species
      * \param interfaceThickness  thickness of the polymer/wall interface
      * \param normalVecId  index of the lattice vector normal to the film
      * \param nMonomer  number of monomer species in the system
      * \param meshDims  grid dimensions of the field in real space
      */
      template <int D>
      __global__ void _hFieldDerivatives(cudaReal ** hDerivs, 
                                         cudaReal const * mask, 
                                         cudaReal const * chiBottom, 
                                         cudaReal const * chiTop, 
                                         cudaReal const interfaceThickness, 
                                         int const normalVecId, 
                                         int const nMonomer, 
                                         dim3 const meshDims)
      {
         // Get position in 3-dimensional grid of threads
         int ix = blockIdx.x * blockDim.x + threadIdx.x;
         int iy = blockIdx.y * blockDim.y + threadIdx.y;
         int iz = blockIdx.z * blockDim.z + threadIdx.z;

         // Get thread index in linearized array
         int id = ix; 
         if (D > 1) id += meshDims.x * iy;
         if (D > 2) id += meshDims.x * meshDims.y * iz;
         
         // Get position d along normalVecId in reduced coordinates
         cudaReal d;
         if (D == 3) {
            if (normalVecId == 2) d = (cudaReal)ix / (cudaReal)meshDims.x;
            if (normalVecId == 1) d = (cudaReal)iy / (cudaReal)meshDims.y;
            if (normalVecId == 0) d = (cudaReal)iz / (cudaReal)meshDims.z;
         } else if (D == 2) {
            if (normalVecId == 1) d = (cudaReal)ix / (cudaReal)meshDims.x;
            if (normalVecId == 0) d = (cudaReal)iy / (cudaReal)meshDims.y;
         } else { // D == 1
            d = (cudaReal)ix / (cudaReal)meshDims.x;
         }

         if (ix < meshDims.x && iy < meshDims.y && iz < meshDims.z) {
            // Get mask value at this gridpoint from global memory
            cudaReal dPhidL = mask[id];

            // Calculate d\phi_m / dL at this gridpoint
            dPhidL = dPhidL * (dPhidL - 1) * 8.0 * (abs(d - 0.5) - 0.5)
                     / interfaceThickness;
            
            // Calculate dh_i / dL at this gridpoint for each species i
            for (int in = 0; in < nMonomer; in++) {
               if (d < 0.5) {
                  hDerivs[in][id] = -1.0 * dPhidL * chiBottom[in];
               } else {
                  hDerivs[in][id] = -1.0 * dPhidL * chiTop[in];
               }
            }
         }
      }

      /*
      * Generate the h fields.
      * 
      * \param hFields  array of h fields (output)
      * \param mask  field containing the mask
      * \param chiBottom  array of chiBottom values for each species
      * \param chiTop  array of chiTop values for each species
      * \param normalVecId  index of the lattice vector normal to the film
      * \param nMonomer  number of monomer species in the system
      * \param meshDims  grid dimensions of the field in real space
      */
      template <int D>
      __global__ void _generateHFields(cudaReal ** hFields, 
                                       cudaReal const * mask, 
                                       cudaReal const * chiBottom, 
                                       cudaReal const * chiTop,
                                       int const normalVecId, 
                                       int const nMonomer, 
                                       dim3 const meshDims)
      {
         // Get position in 3-dimensional grid of threads
         int ix = blockIdx.x * blockDim.x + threadIdx.x;
         int iy = blockIdx.y * blockDim.y + threadIdx.y;
         int iz = blockIdx.z * blockDim.z + threadIdx.z;

         // Get thread index in linearized array
         int id = ix; 
         if (D > 1) id += meshDims.x * iy;
         if (D > 2) id += meshDims.x * meshDims.y * iz;

         // Determine if this grid point is in the "top" half
         bool topHalf;
         if (D == 3) {
            if (normalVecId == 2) topHalf = (ix > meshDims.x / 2);
            if (normalVecId == 1) topHalf = (iy > meshDims.y / 2);
            if (normalVecId == 0) topHalf = (iz > meshDims.z / 2);
         } else if (D == 2) {
            if (normalVecId == 1) topHalf = (ix > meshDims.x / 2);
            if (normalVecId == 0) topHalf = (iy > meshDims.y / 2);
         } else { // D == 1
            topHalf = (ix > meshDims.x / 2);
         }

         if (ix < meshDims.x && iy < meshDims.y && iz < meshDims.z) {
            // Get mask value at this gridpoint from global memory
            cudaReal maskVal = mask[id];
            
            // Calculate h field value at this gridpoint for each species i
            for (int in = 0; in < nMonomer; in++) {
               if (topHalf) {
                  hFields[in][id] = (1.0 - maskVal) * chiTop[in];
               } else {
                  hFields[in][id] = (1.0 - maskVal) * chiBottom[in];
               }
            }
         }
      }

   }

   /*
   * Default constructor
   */
   template <int D>
   ExtGenFilm<D>::ExtGenFilm()
    : ExtGenFilmBase<D>::ExtGenFilmBase(),
      sysPtr_(0)
   {  setClassName("ExtGenFilm"); }
   
   /*
   * Constructor
   */
   template <int D>
   ExtGenFilm<D>::ExtGenFilm(System<D>& sys)
    : ExtGenFilmBase<D>::ExtGenFilmBase(),
      sysPtr_(&sys)
   {  setClassName("ExtGenFilm"); }

   /*
   * Destructor
   */
   template <int D>
   ExtGenFilm<D>::~ExtGenFilm()
   {}

   /*
   * Get contribution to the stress from these external fields
   */
   template <int D>
   double ExtGenFilm<D>::stressTerm(int paramId) const
   {
      // If walls are athermal then there is no external field, so no
      // contribution to the stress.
      if (isAthermal()) return 0.0;
      
      // If paramId is not normalVecId, there is no stress contribution
      UTIL_CHECK(normalVecId() >= 0); // normalVecId has been set
      int nvParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                                 system().domain().lattice());
      if (nvParamId != paramId) return 0.0;

      // If this point is reached, calculate the stress contribution
      // from the external fields.
      UTIL_CHECK(isGenerated());
      UTIL_CHECK(interfaceThickness() > 0); 
      UTIL_CHECK(system().hasMask());
      UTIL_CHECK(system().hasExternalFields());

      // Setup
      int nMonomer = system().mixture().nMonomer();

      DArray< RField<D> > hDerivatives;
      HostDArray<cudaReal*> hDerivPtrs_h(nMonomer);
      DeviceArray<cudaReal*> hDerivPtrs(nMonomer);

      hDerivatives.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         hDerivatives[i].allocate(system().mesh().dimensions());
         hDerivPtrs_h[i] = hDerivatives[i].cArray();
      }
      hDerivPtrs = hDerivPtrs_h;

      HostDArray<cudaReal>  chiBottom_h(nMonomer), chiTop_h(nMonomer);
      DeviceArray<cudaReal> chiBottom_d(nMonomer), chiTop_d(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         chiBottom_h[i] = chiBottom(i);
         chiTop_h[i] = chiTop(i);
      }
      chiBottom_d = chiBottom_h;
      chiTop_d = chiTop_h;

      // Set D-dimensional GPU configuration
      ThreadMesh::setConfig(system().mesh().dimensions(), true);
      dim3 gridDims = ThreadMesh::gridDims();
      dim3 blockDims = ThreadMesh::blockDims();

      // Compute the derivative of h fields with respect to film thickness
      _hFieldDerivatives<D><<<gridDims, blockDims>>>
         (hDerivPtrs.cArray(), system().mask().rgrid().cArray(),
          chiBottom_d.cArray(), chiTop_d.cArray(), interfaceThickness(),
          normalVecId(), nMonomer, ThreadMesh::meshDims());
      
      // Integrate resulting fields to get the stress
      double stress;
      for (int i = 0; i < nMonomer; i++) {
         stress += Reduce::innerProduct(hDerivatives[i], system().c().rgrid(i));
      }
      stress /= system().mask().phiTot() * hDerivatives[0].capacity();
      return stress;
   }

   /*
   * Allocate container necessary to generate and store field
   */
   template <int D>
   void ExtGenFilm<D>::allocate()
   {
      UTIL_CHECK(system().unitCell().isInitialized());

      // Make sure h field container has access to a fieldIo
      system().h().setFieldIo(system().domain().fieldIo());

      // Allocate the external field containers if needed
      if (!isAthermal()) {
         if (!system().h().isAllocatedRGrid()) {
            system().h().allocateRGrid(system().mesh().dimensions());
         }
         if (system().iterator().isSymmetric()) {
            UTIL_CHECK(system().basis().isInitialized());
            if (!system().h().isAllocatedBasis()) {
               system().h().allocateBasis(system().basis().nBasis());
            }
         }
      }
   }

   /**
   * Generate the fields and store where the Iterator can access.
   */
   template <int D>
   void ExtGenFilm<D>::generate()
   {
      // If walls are athermal then there is no external field needed.
      // If an external field already exists in the System, we need to
      // overwrite it with a field of all zeros, otherwise do nothing
      if ((isAthermal()) && (!isGenerated())) return;

      // If this point is reached, external field must be generated
      UTIL_CHECK(system().h().isAllocatedRGrid());
      UTIL_CHECK(system().mask().isAllocatedRGrid());
      if (system().iterator().isSymmetric()) {
         UTIL_CHECK(system().h().isAllocatedBasis());
         UTIL_CHECK(system().mask().isAllocatedBasis());
      }
      
      UTIL_CHECK(normalVecId() >= 0);

      // Set chiBottomCurrent_, chiTopCurrent_, and parametersCurrent_
      chiBottomCurrent_ = chiBottom();
      chiTopCurrent_ = chiTop();
      normalVecCurrent_ = systemLatticeVector(normalVecId());

      // Setup
      int nMonomer = system().mixture().nMonomer();

      DArray< RField<D> > hFields;
      HostDArray<cudaReal*> hPtrs_h(nMonomer);
      DeviceArray<cudaReal*> hPtrs(nMonomer);

      hFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         hFields[i].allocate(system().mesh().dimensions());
         hPtrs_h[i] = hFields[i].cArray();
      }
      hPtrs = hPtrs_h;

      HostDArray<cudaReal>  chiBottom_h(nMonomer), chiTop_h(nMonomer);
      DeviceArray<cudaReal> chiBottom_d(nMonomer), chiTop_d(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         chiBottom_h[i] = chiBottom(i);
         chiTop_h[i] = chiTop(i);
      }
      chiBottom_d = chiBottom_h;
      chiTop_d = chiTop_h;

      // Set D-dimensional GPU configuration
      ThreadMesh::setConfig(system().mesh().dimensions(), true);
      dim3 gridDims = ThreadMesh::gridDims();
      dim3 blockDims = ThreadMesh::blockDims();

      // Generate the h fields
      _generateHFields<D><<<gridDims, blockDims>>>
         (hPtrs.cArray(), system().mask().rgrid().cArray(), 
          chiBottom_d.cArray(), chiTop_d.cArray(), normalVecId(), 
          nMonomer, ThreadMesh::meshDims());

      // Pass h into the System
      system().h().setRGrid(hFields, system().iterator().isSymmetric());
   }

}
}

#endif
