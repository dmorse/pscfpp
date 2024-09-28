#ifndef RPC_EXT_GEN_FILM_TPP
#define RPC_EXT_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExtGenFilm.h"
#include <rpc/field/FieldIo.h>
#include <prdc/cpu/RField.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FArray.h>

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
   * Allocate container necessary to generate and store field
   */
   template <int D>
   void ExtGenFilm<D>::allocate()
   {
      UTIL_CHECK(system().domain().basis().isInitialized());
      UTIL_CHECK(system().domain().unitCell().isInitialized());

      system().h().setFieldIo(system().domain().fieldIo());

      // Allocate the external field containers if needed
      if (!isAthermal()) {
         if (!system().h().isAllocatedRGrid()) {
            system().h().allocateRGrid(system().domain().mesh().dimensions());
         }
         if (!system().h().isAllocatedBasis()) {
            system().h().allocateBasis(system().domain().basis().nBasis());
         }
      }
   }

   /**
   * Generate the fields and store where the Iterator can access.
   */
   template <int D>
   void ExtGenFilm<D>::generate()
   {
      // Make sure normalVecId_ is set
      if (normalVecId_ < 0) maskNormalVecId();

      // Set chiBottomCurrent_, chiTopCurrent_, and parametersCurrent_
      chiBottomCurrent_ = chiBottom();
      chiTopCurrent_ = chiTop();
      parametersCurrent_ = systemLatticeParameters();

      // If walls are athermal then there is no external field needed.
      // If an external field already exists in the System, we need to
      // overwrite it with a field of all zeros, otherwise do nothing
      if ((isAthermal()) && (!isGenerated())) return;

      // If this point is reached, external field must be generated
      UTIL_CHECK(system().h().isAllocatedRGrid());
      UTIL_CHECK(system().h().isAllocatedBasis());
      UTIL_CHECK(system().mask().isAllocated());

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
                  if (coords[normalVecId_] < (dim[normalVecId_]/2)) {
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
      system().h().setRGrid(hRGrid,true);
   }

   /**
   * Use the mask to determine the value of normalVecId
   */
   template <int D>
   void ExtGenFilm<D>::maskNormalVecId()
   {
      // If normalVecId_ has already been set, do nothing.
      if (normalVecId_ >= 0) return;
      
      // Determine normalVecId by checking the value of the mask in the 
      // middle of each edge of the unit cell (one edge per dimension). 
      // Dimension normalVecId will have a mask value around 1, the others
      // will have a mask value around 0.

      double maxVal = 0;

      // Make IntVec position array that denotes the origin (D zeros)
      IntVec<D> position;
      for (int i = 0; i < D; i++)
         position[i] = 0;
      
      // pointer to mask RField
      RField<D> const & maskPtr = system().mask().rgrid(); 

      for (int i = 0; i < D; i++) {
         position[i] = system().domain().mesh().dimension(i) / 2;
         int rank = system().domain().mesh().rank(position);
         if (maskPtr[rank] > maxVal) {
            maxVal = maskPtr[rank];
            normalVecId_ = i;
         }
         position[i] = 0;
      }

      // Make sure value was actually set and maxVal makes sense
      UTIL_CHECK(normalVecId_ >= 0);
      UTIL_CHECK(maxVal > 0.99);
   } 
}
}

#endif
