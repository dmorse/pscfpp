#ifndef RPC_EXT_GEN_FILM_TPP
#define RPC_EXT_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExtGenFilm.h"
#include "Iterator.h"
#include <rpc/field/FieldIo.h>
#include <prdc/cpu/RField.h>
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
      int nx = system().domain().mesh().size();
      RField<D> const & maskRGrid = system().mask().rgrid();
      RField<D> maskDeriv, hDeriv;
      maskDeriv.allocate(system().domain().mesh().dimensions());
      hDeriv.allocate(system().domain().mesh().dimensions());
      FArray<int,3> coords;
      int x, y, z;
      int counter = 0;
      double d, maskVal;
      double term = 0.0;

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

      // Compute the derivative of the mask with respect to film thickness
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

               maskDeriv[counter] = maskVal * (maskVal - 1) * 8.0
                                    * (std::abs(d - 0.5) - 0.5)
                                    / interfaceThickness();
               counter++;
            }
         }
      }

      for (int i = 0; i < nMonomer; i++) {
         // Compute the h field derivative with respect to film thickness
         counter = 0;
         for (x = 0; x < dim[0]; x++) {
            coords[0] = x;
            for (y = 0; y < dim[1]; y++) {
               coords[1] = y;
               for (z = 0; z < dim[2]; z++) {
                  coords[2] = z;
                  if (coords[normalVecId()] < (dim[normalVecId()] / 2)) {
                     hDeriv[counter] = -1.0 * maskDeriv[counter] * chiBottom(i);
                  } else {
                     hDeriv[counter] = -1.0 * maskDeriv[counter] * chiTop(i);
                  }
                  counter++;
               }
            }
         }

         // Get the integral term in the stress
         RField<D> const & c = system().c().rgrid(i);
         for (int i = 0; i < nx; i++) {
            term += c[i] * hDeriv[i];
         }
      }
      term /= (system().mask().phiTot() * nx);
      return term;
   }

   /*
   * Allocate container necessary to generate and store field
   */
   template <int D>
   void ExtGenFilm<D>::allocate()
   {
      UTIL_CHECK(system().domain().unitCell().isInitialized());

      // Make sure h field container has access to a fieldIo
      system().h().setFieldIo(system().domain().fieldIo());

      // Allocate the external field containers if needed
      if (!isAthermal()) {
         if (!system().h().isAllocatedRGrid()) {
            system().h().allocateRGrid(system().domain().mesh().dimensions());
         }
         if (system().iterator().isSymmetric()) {
            UTIL_CHECK(system().domain().basis().isInitialized());
            if (!system().h().isAllocatedBasis()) {
               system().h().allocateBasis(system().domain().basis().nBasis());
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
      // Set chiBottomCurrent_, chiTopCurrent_, and parametersCurrent_
      chiBottomCurrent_ = chiBottom();
      chiTopCurrent_ = chiTop();
      normalVecCurrent_ = systemLatticeVector(normalVecId());

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

      int nm = systemNMonomer();

      // Create a 3 element vector 'dim' that contains the grid dimensions.
      // If system is 2D (1D), then the z (y & z) dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimension(ind);
         } else {
            dim[ind] = 1;
         }
      }

      // Get pointer to mask RField
      RField<D> const & maskPtr = system().mask().rgrid();

      // Generate an r-grid representation of the external fields
      DArray< RField<D> > hRGrid;
      hRGrid.allocate(nm);
      for (int i = 0; i < nm; i++) {
         hRGrid[i].allocate(system().domain().mesh().dimensions());
      }

      int i, x, y, z;
      int counter = 0;
      FArray<int,3> coords;
      double rhoW;

      for (i = 0; i < nm; i++) {
         for (x = 0; x < dim[0]; x++) {
            coords[0] = x;
            for (y = 0; y < dim[1]; y++) {
               coords[1] = y;
               for (z = 0; z < dim[2]; z++) {
                  coords[2] = z;

                  // Calculate wall volume fraction (rho_w) at gridpoint 
                  // (x,y,z)
                  rhoW = maskPtr[counter];
                  if (coords[normalVecId()] < (dim[normalVecId()]/2)) {
                     hRGrid[i][counter++] = (1.0-rhoW) * chiBottom(i);
                  } else {
                     hRGrid[i][counter++] = (1.0-rhoW) * chiTop(i);
                  }
               }
            }
         }
         counter = 0;
      } 

      // Pass h into the System
      system().h().setRGrid(hRGrid, system().iterator().isSymmetric());
   }

}
}

#endif
