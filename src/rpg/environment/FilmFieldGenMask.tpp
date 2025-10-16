#ifndef RPG_MASK_GEN_FILM_TPP
#define RPG_MASK_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenMask.h"
#include <rpg/scft/ScftThermo.h>
#include <rpg/scft/iterator/Iterator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/paramIdConversions.h>
#include <prdc/cuda/resources.h>
#include <pscf/inter/Interaction.h>
#include <cmath>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   // CUDA kernels (anonymous namespace - used only in this file)
   namespace {

      /*
      * Compute the derivative of mask field with respect to film thickness.
      *
      * (note: indices are inverted from their RField representation. So,
      * for a 3D system, x, y, and z in an RField become z, y, and x in the 
      * CUDA thread grid, respectively. In 2D, x, y becomes y, x. This 
      * allows for memory coalescing; see ThreadMesh::setConfig 
      * documentation for more information.)
      *
      * \param deriv  field containing the mask derivative (output)
      * \param mask  field containing the mask
      * \param interfaceThickness  thickness of the polymer/wall interface
      * \param normalVecId  index of the lattice vector normal to the film
      * \param meshDims  grid dimensions of the field in real space
      */
      template <int D>
      __global__ void _maskDerivative(cudaReal * deriv, 
                                      cudaReal const * mask, 
                                      cudaReal const interfaceThickness, 
                                      int const normalVecId, 
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
            cudaReal maskVal = mask[id];

            // Calculate d\phi_m / dL at this gridpoint
            deriv[id] = maskVal * (maskVal - 1) * 8.0 * (abs(d - 0.5) - 0.5)
                        / interfaceThickness;
         }
      }

      /*
      * Generate the mask.
      *
      * \param mask  field containing the mask (output)
      * \param nvLength  length of the lattice vector normal to the film
      * \param excludedThickness  thickness of region occupied by the wall
      * \param interfaceThickness  thickness of the polymer/wall interface
      * \param normalVecId  index of the lattice vector normal to the film
      * \param meshDims  grid dimensions of the field in real space
      */
      template <int D>
      __global__ void _generateMask(cudaReal * mask, cudaReal const nvLength,
                                    cudaReal const interfaceThickness,
                                    cudaReal const excludedThickness,
                                    int const normalVecId, dim3 const meshDims)
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

         // Convert d to unreduced coordinates
         d *= nvLength;

         if (ix < meshDims.x && iy < meshDims.y && iz < meshDims.z) {
            // Calculate wall volume fraction at this gridpoint
               mask[id] =  
                  (0.5 * (1.0 - tanh(
                     (
                        0.5 * (excludedThickness - nvLength) +
                        abs(d - (nvLength / 2.0))
                     ) 
                     * 4.0 / interfaceThickness
                  )));
         }
      }

   }

   /*
   * Default constructor
   */
   template <int D>
   FilmFieldGenMask<D>::FilmFieldGenMask()
    : FilmFieldGenMaskBase<D>::FilmFieldGenMaskBase(),
      sysPtr_(0)
   {  setClassName("FilmFieldGenMask"); }
   
   /*
   * Constructor
   */
   template <int D>
   FilmFieldGenMask<D>::FilmFieldGenMask(System<D>& sys)
    : FilmFieldGenMaskBase<D>::FilmFieldGenMaskBase(),
      sysPtr_(&sys)
   {  setClassName("FilmFieldGenMask"); }

   /*
   * Destructor
   */
   template <int D>
   FilmFieldGenMask<D>::~FilmFieldGenMask()
   {}

   /*
   * Get contribution to the stress from this mask
   */
   template <int D>
   double FilmFieldGenMask<D>::stress(int paramId) const
   {
      // If paramId is not normalVecId, there is no stress contribution
      UTIL_CHECK(normalVecId() >= 0);
      int normalVecParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                                system().domain().lattice());
      if (normalVecParamId != paramId) return 0.0;

      // If this point is reached, stress contribution must be calculated
      UTIL_CHECK(system().mask().hasData());
      
      // Get the length of the lattice basis vector normal to the walls
      double nvLength = system().domain().unitCell().parameter(paramId);

      // Get the volume fraction of the unit cell occupied by polymers
      double phiTot = system().mask().phiTot();

      // Create the derivative field in rgrid format
      RField<D> deriv;
      deriv.allocate(system().domain().mesh().dimensions());
      
      // Set D-dimensional GPU configuration
      ThreadMesh::setConfig(system().domain().mesh().dimensions(), true);
      dim3 gridDims = ThreadMesh::gridDims();
      dim3 blockDims = ThreadMesh::blockDims();

      // Compute derivative of the mask with respect to film thickness
      _maskDerivative<D><<<gridDims, blockDims>>>
         (deriv.cArray(), system().mask().rgrid().cArray(),
          interfaceThickness(), normalVecId(), ThreadMesh::meshDims());
      
      // Get xi, the Lagrange multiplier field, in rgrid format
      RField<D> xi;
      xi.allocate(system().domain().mesh().dimensions());

      if (system().h().hasData()) {
         VecOp::subVV(xi, system().w().rgrid(0), system().h().rgrid(0));
      } else {
         VecOp::eqV(xi, system().w().rgrid(0));
      }

      int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; i++) {
         double chi = system().interaction().chi(0, i);
         if (fabs(chi) > 1e-6) { // if chi is nonzero
            VecOp::addEqVc(xi, system().c().rgrid(i), -1.0 * chi);
         }
      }

      // Integrate deriv * xi to get the integral term in the stress
      double intTerm = Reduce::innerProduct(deriv, xi);
      intTerm /= phiTot * deriv.capacity();

      // Get the pressure term in the stress
      if (!sysPtr_->scft().hasData()) {
         sysPtr_->scft().compute();
      }
      double pSys = sysPtr_->scft().pressure();
      double pTerm = pSys * excludedThickness() / 
                     (phiTot * nvLength * nvLength);
      
      return pTerm - intTerm;
   }

   template <int D>
   double FilmFieldGenMask<D>::modifyStress(int paramId, double stress) 
   const
   {
      int nvParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                             system().domain().lattice());

      // If paramId is not normalVecId, there is no stress modification
      if (nvParamId != paramId) return stress;

      // If this point is reached, stress must be modified
      UTIL_CHECK(system().mask().hasData());

      if (!hasFBulk()) {
         UTIL_THROW("fBulk must be set before calling modifyStress.");
      }

      // Get system free energy
      if (!sysPtr_->scft().hasData()) {
         sysPtr_->scft().compute();
      }
      double fSys = sysPtr_->scft().fHelmholtz();

      // Get the length L of the basis vector normal to the walls
      double L = system().domain().unitCell().parameter(paramId);

      double modifiedStress = (stress * (L - excludedThickness())) 
                              + (fSys - fBulk_);

      return modifiedStress;
   }

   /*
   * Compute the field and store where the System can access
   */
   template <int D>
   void FilmFieldGenMask<D>::compute()
   {
      UTIL_CHECK(normalVecId() >= 0);
      UTIL_CHECK(interfaceThickness() > 0);
      UTIL_CHECK(excludedThickness() > interfaceThickness());
      if (system().iterator().isSymmetric()) {
         UTIL_CHECK(system().domain().basis().isInitialized());
      }

      // Get the length of the lattice basis vector normal to the walls
      int paramId;
      paramId = convertFullParamIdToReduced<D>(normalVecId(),
                                               system().domain().lattice());
      double nvLength = system().domain().unitCell().parameter(paramId);

      // Setup
      RField<D> rGrid;
      rGrid.allocate(system().domain().mesh().dimensions());

      // Set D-dimensional GPU configuration
      ThreadMesh::setConfig(system().domain().mesh().dimensions(), true);
      dim3 gridDims = ThreadMesh::gridDims();
      dim3 blockDims = ThreadMesh::blockDims();

      // Generate the mask
      _generateMask<D><<<gridDims, blockDims>>>
         (rGrid.cArray(), nvLength, interfaceThickness(), 
          excludedThickness(), normalVecId(), ThreadMesh::meshDims());

      // Store this mask in System
      system().mask().setRGrid(rGrid,system().iterator().isSymmetric());

      // Store lattice vector normal to film used to construct this mask
      normalVecCurrent_ = systemLatticeVector(normalVecId());
   }

   /*
   * Sets flexible lattice parameters to be compatible with the mask.
   */
   template <int D>
   void FilmFieldGenMask<D>::setFlexibleParams() const
   {
      if (system().iterator().isFlexible()) {
         FSArray<bool,6> updated;
         updated = modifyFlexibleParams(system().iterator().flexibleParams(),
                                        system().domain().unitCell());
         sysPtr_->iterator().setFlexibleParams(updated);
      }
   }

} // namespace Rpg
} // namespace Pscf
#endif
