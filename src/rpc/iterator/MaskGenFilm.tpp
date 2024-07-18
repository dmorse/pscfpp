#ifndef RPC_MASK_GEN_FILM_TPP
#define RPC_MASK_GEN_FILM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilm.h"
#include <prdc/cpu/RField.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/math/RealVec.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FArray.h>
#include <rpc/field/FieldIo.h>

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
   MaskGenFilm<D>::MaskGenFilm()
    : MaskGenFilmBase<D>::MaskGenFilmBase(),
      sysPtr_(0)
   {  setClassName("MaskGenFilm"); }
   
   /*
   * Constructor
   */
   template <int D>
   MaskGenFilm<D>::MaskGenFilm(System<D>& sys)
    : MaskGenFilmBase<D>::MaskGenFilmBase(),
      sysPtr_(&sys)
   {  setClassName("MaskGenFilm"); }

   /*
   * Destructor
   */
   template <int D>
   MaskGenFilm<D>::~MaskGenFilm()
   {}

   /*
   * Allocate container necessary to generate and store field
   */
   template <int D>
   void MaskGenFilm<D>::allocate()
   {
      UTIL_CHECK(system().basis().isInitialized());
      UTIL_CHECK(system().unitCell().isInitialized());

      system().mask().setFieldIo(system().fieldIo());

      // Allocate the mask containers if needed
      if (!system().mask().isAllocated()) {
         system().mask().allocate(system().basis().nBasis(), 
                                  system().mesh().dimensions());
      }
   }

   /*
   * Generate the field and store where the Iterator can access
   */
   template <int D>
   void MaskGenFilm<D>::generate()
   {
      UTIL_CHECK(interfaceThickness() > 0);
      UTIL_CHECK(excludedThickness() > interfaceThickness());
      UTIL_CHECK(system().mask().isAllocated());

      // Get the length L of the lattice basis vector normal to the walls
      RealVec<D> a;
      a = system().domain().unitCell().rBasis(normalVecId());
      double normSqd(0.0); // norm squared
      for (int i = 0; i < D; i++) {
         normSqd += a[i]*a[i];
      }
      double L(sqrt(normSqd));

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
      system().mask().setRGrid(rGrid,true);

      // Store lattice parameters associated with this maskBasis
      parameters_ = system().domain().unitCell().parameters();

   }

   // Explicit Specializations for setFlexibleParams are in MaskGenFilm.cpp
}
}
#endif