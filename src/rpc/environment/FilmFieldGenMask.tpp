#ifndef RPC_MASK_GEN_FILM_TPP
#define RPC_MASK_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenMask.h"
#include <rpc/field/FieldIo.h>
#include <rpc/scft/iterator/Iterator.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/paramIdConversions.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FArray.h>
#include <cmath>

namespace Pscf {
namespace Rpc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Default constructor
   */
   template <int D>
   FilmFieldGenMask<D>::FilmFieldGenMask()
    : FilmFieldGenMaskBase<D>::FilmFieldGenMaskBase(),
      sysPtr_(nullptr)
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
      UTIL_CHECK(sysPtr_);
      int normalVecParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                                system().domain().lattice());

      // If paramId is not normalVecId, there is no stress contribution
      if (normalVecParamId != paramId) return 0.0;

      // If this point is reached, stress contribution must be calculated
      UTIL_CHECK(system().hasMask());
      
      // Get the length of the lattice basis vector normal to the walls
      double nvLength = system().domain().unitCell().parameter(paramId);

      // Get the volume fraction of the unit cell occupied by polymers
      double phiTot = system().mask().phiTot();

      // Create a 3 element vector 'dim' that contains the grid dimensions.
      // If system is 2D (1D), then the z (y & z) dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimensions()[ind];
         } else {
            dim[ind] = 1;
         }
      }

      // Create the derivative field in rgrid format
      RField<D> deriv;
      FArray<int,3> coords;
      int x, y, z;
      int counter = 0;
      double d, maskVal;

      deriv.allocate(system().domain().mesh().dimensions());
      RField<D> const & maskRGrid = system().mask().rgrid();
      
      for (x = 0; x < dim[0]; x++) {
         coords[0] = x;
         for (y = 0; y < dim[1]; y++) {
            coords[1] = y;
            for (z = 0; z < dim[2]; z++) {
               coords[2] = z;
               maskVal = maskRGrid[counter];

               // Get the distance 'd' traveled along the lattice basis 
               // vector normal to the walls, in reduced coordinates
               d = (double)coords[normalVecId()] / 
                     (double)dim[normalVecId()];

               deriv[counter] = maskVal * (maskVal - 1) * 8.0
                                 * (std::abs(d - 0.5) - 0.5)
                                 / interfaceThickness();
               counter++;
            }
         }
      }

      // Get xi, the Lagrange multiplier field, in rgrid format
      int nMonomer = system().mixture().nMonomer();
      int nx = system().domain().mesh().size();
      double chi;
      RField<D> xi;
      xi.allocate(system().domain().mesh().dimensions());

      for (int i = 0; i < nx; i++) {
         xi[i] = system().w().rgrid(0)[i];
      }
      
      for (int in = 0; in < nMonomer; in++) {
         chi = system().interaction().chi(0,in);
         if (fabs(chi) > 1e-6) { // if chi is nonzero
            for (int i = 0; i < nx; i++) {
               xi[i] -= system().c().rgrid(in)[i] * chi;
            }
         }
      }

      if (system().hasExternalFields()) {
         for (int i = 0; i < nx; i++) {
            xi[i] -= system().h().rgrid(0)[i];
         }
      }

      // Get the integral term in the stress
      double intTerm = 0.0;
      for (int i = 0; i < nx; i++) {
         intTerm += xi[i] * deriv[i];
      }
      intTerm /= (phiTot * nx);

      // Get the pressure term in the stress
      if (!sysPtr_->hasFreeEnergy()) {
         sysPtr_->computeFreeEnergy();
      }
      double pSys = sysPtr_->pressure();
      double pTerm = pSys * excludedThickness() / 
                     (phiTot * nvLength * nvLength);
      
      return pTerm - intTerm;
   }

   template <int D>
   double FilmFieldGenMask<D>::modifyStress(int paramId, double stress) 
   const
   {
      UTIL_CHECK(sysPtr_);
      int nvParamId = convertFullParamIdToReduced<D>(normalVecId(),
                                             system().domain().lattice());
      if (nvParamId == paramId) {
         UTIL_CHECK(system().hasMask());

         if (!hasFBulk()) {
            UTIL_THROW("fBulk must be set before calculating this stress.");
         }

         // Get system free energy
         if (!sysPtr_->hasFreeEnergy()) {
            sysPtr_->computeFreeEnergy();
         }
         double fSys = sysPtr_->fHelmholtz();

         // Get the length L of the basis vector normal to the walls
         double L = system().domain().unitCell().parameter(paramId);

         double modifiedStress = (stress * (L - excludedThickness())) 
                                 + (fSys - fBulk_);

         return modifiedStress;
      } else {
         return stress;
      }
   }

   /*
   * Compute the field and store where the System can access
   */
   template <int D>
   void FilmFieldGenMask<D>::compute()
   {
      UTIL_CHECK(sysPtr_);
      UTIL_CHECK(interfaceThickness() > 0);
      UTIL_CHECK(excludedThickness() > interfaceThickness());
      UTIL_CHECK(normalVecId() >= 0);
      if (system().iterator().isSymmetric()) {
         UTIL_CHECK(system().domain().basis().isInitialized());
      }
      
      // Get the length L of the lattice basis vector normal to the walls
      int paramId = convertFullParamIdToReduced<D>(normalVecId(),
                                                system().domain().lattice());
      double L = system().domain().unitCell().parameter(paramId);

      // Create a 3 element vector 'dim' that contains the grid dimensions.
      // If system is 2D (1D), then the z (y & z) dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimensions()[ind];
         } else {
            dim[ind] = 1;
         }
      }

      // Generate an r-grid representation of the walls
      RField<D> rGrid;
      rGrid.allocate(system().domain().mesh().dimensions());
      int x, y, z;
      int counter = 0;
      FArray<int,3> coords;
      double d, rhoW;

      for (x = 0; x < dim[0]; x++) {
         coords[0] = x;
         for (y = 0; y < dim[1]; y++) {
            coords[1] = y;
            for (z = 0; z < dim[2]; z++) {
               coords[2] = z;

               // Get the distance 'd' traveled along the lattice basis 
               // vector that is normal to the walls
               d = coords[normalVecId()] * L / dim[normalVecId()];

               // Calculate wall volume fraction (rhoW) at gridpoint (x,y,z)
               rhoW = 0.5 * (1 + tanh(4 * (((.5 * (excludedThickness()-L)) + 
                                  fabs(d - (L/2))) / interfaceThickness())));
               rGrid[counter++] = 1-rhoW;
            }
         }
      }

      // Store this mask in System
      system().mask().setRGrid(rGrid, system().iterator().isSymmetric());

      // Store lattice vector normal to film used to construct this mask
      normalVecCurrent_ = systemLatticeVector(normalVecId());

   }

   /*
   * Sets flexible lattice parameters to be compatible with the mask.
   */
   template <int D>
   void FilmFieldGenMask<D>::setFlexibleParams() const
   {
      UTIL_CHECK(sysPtr_);
      if (system().iterator().isFlexible()) {
         FSArray<bool,6> updated;
         updated = modifyFlexibleParams(system().iterator().flexibleParams(),
                                        system().domain().unitCell());
         sysPtr_->iterator().setFlexibleParams(updated);
      }
   }
   
} // namespace Rpc
} // namespace Pscf
#endif
