#ifndef PSSP_R_MESH_FIELD_H
#define PSSP_R_MESH_FIELD_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pssp
{

   using namespace Util;
   using namespace Pscf;

   /**
   * Field of real double precision values on an FFT mesh.
   */
   class RMeshField : public Field<double>
   {

   public:

      /**
      * Default constructor.
      */
      RMeshField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RMeshField to be copied.
      */
      RMeshField(const RMeshField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RMeshField();

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RMeshField is already allocated.
      *
      * \param dimensions vector containing number of grid points in each direction.
      */
      template <int D>
      void allocate(const IntVec<D>& meshDimensions);

      /**
      * Get the dimensions of the grid for which this was allocated.
      *
      * \throw Exception if dimensions of space do not match.
      *
      * \param dimensions vector containing number of grid points in each direction.
      */
      template <int D>
      void getMeshDimension(IntVec<D>& meshDimensions);

   private:

      // Dimension of space (1, 2, or 3)
      int spaceDimension_;

      // Vector containing number of grid points in each direction.
      IntVec<3> meshDimensions_;

   };

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RMeshField::allocate(const IntVec<D>& meshDimensions)
   {
      // Preconditions
      UTIL_CHECK(D > 0);
      UTIL_CHECK(D < 4);

      // Initialize mesh dimensions to zero
      for (int i = 0; i < 3; ++i) {
         meshDimensions_[i] = 1;
      }

      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         size *= meshDimensions[i];
         meshDimensions_[i] = meshDimensions[i];
      }
      spaceDimension_ = D;
      Field<double>::allocate(size);
   }

   /**
   * Get the dimensions of the grid for which this was allocated.
   */
   template <int D>
   void RMeshField::getMeshDimension(IntVec<D>& meshDimensions)
   {
      if (D != spaceDimension_) {
         UTIL_THROW("Argument with wrong number of spatial dimensions");
      } else {
         for (int i = 0; i < D; ++i) {
            meshDimensions[i] = meshDimensions_[i];
         }
      }
   }

}
#endif
