#ifndef PRDC_MASK_GEN_FILM_BASE_H
#define PRDC_MASK_GEN_FILM_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/FieldGenerator.h>  // Base class
#include <pscf/math/RealVec.h>        // container
#include <util/containers/FSArray.h>  // container
#include <iostream>
#include <string>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base class defining mask that imposes thin film confinement.
   * 
   * This is a base class for MaskGenFilm that defines all traits of a 
   * MaskGenFilm that do not require access to the System (System access is
   * needed, for example, to get the space group and set the mask field).
   * 
   * If the user chooses a MaskGenFilm object to construct the mask, then 
   * the system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell. The
   * shape of the mask is defined by three input parameters: normalVecId,
   * excludedThickness, and interfaceThickness. See \ref 
   * user_thin_films_page for more information. 
   * 
   * \ingroup Prdc_Iterator_Module
   */
   template <int D>
   class MaskGenFilmBase : public FieldGenerator
   {

   public:

      /**
      * Constructor
      * 
      * \param itr  Iterator parent object
      */
      MaskGenFilmBase();

      /**
      * Destructor
      */
      ~MaskGenFilmBase();

      /**
      * Read and initialize.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Check that the system is compatible with this field.
      * 
      * This method calls setFlexibleParams, checkLatticeVectors, and
      * checkSpaceGroup.
      */
      void checkCompatibility();

      /**
      * Check whether system has changed such that the field needs updating.
      */
      bool updateNeeded() const;

      /**
      * Get value of normalVecId.
      */
      int normalVecId() const;

      /**
      * Get value of interfaceThickness.
      */
      double interfaceThickness() const;

      /**
      * Get value of excludedThickness.
      */
      double excludedThickness() const;

      /**
      * Check whether the field has been generated.
      */
      bool isGenerated() const = 0;

   protected:

      /**
      * Check that space group is compatible with the mask.
      */
      void checkSpaceGroup() const;

      /**
      * Check that lattice vectors are compatible with thin film constraint.
      * 
      * Check that user-defined lattice basis vectors (stored in the
      * Domain member of the parent System object) are compatible with 
      * thin film confinement. The lattice basis vector with index 
      * normalVecId should be normal to the walls, while any other lattice
      * basis vectors must be parallel to the walls.
      */
      virtual void checkLatticeVectors() const;

      /**
      * Allocate container necessary to generate and store field.
      */ 
      virtual void allocate() = 0;

      /**
      * Generate the field and store where the Iterator can access.
      */
      virtual void generate() = 0;

      /**
      * Sets flexible lattice parameters to be compatible with the mask.
      * 
      * An iterator for a thin film SCFT calculation should allow for
      * some lattice parameters to be fixed, while others are held 
      * constant. Subclasses should define this method to set these
      * lattice parameters appropriately so that the film thickness is
      * held constant, while the other lattice parameters are allowed
      * to be flexible if the user chooses.
      */
      virtual void setFlexibleParams() = 0;

      /**
      * Get the space group name for this system.
      */
      virtual std::string systemSpaceGroup() const = 0;

      /**
      * Get the lattice parameters for this system.
      */
      virtual FSArray<double, 6> systemLatticeParameters() const = 0;

      /**
      * Get one of the lattice vectors for this system.
      * 
      * \param id  index of the desired lattice parameter
      */
      virtual RealVec<D> systemLatticeVector(int id) const = 0;

      /**
      * The lattice parameters used to generate the current external fields.
      * 
      * This array is set to be equal to the system lattice parameters each
      * time the external fields are generated. The system's lattice 
      * parameters may then change, and this parametersCurrent_ array is
      * used to detect whether they have changed. This is used to determine 
      * whether a new set of external fields needs to be generated.
      */
      FSArray<double, 6> parametersCurrent_;

   private:

      /// Lattice basis vector that is normal to the walls
      int normalVecId_;

      /// Interface thickness
      double interfaceThickness_;

      /// Excluded (wall) thickness
      double excludedThickness_;

      using FieldGenerator::type_;

   };

   // Inline member functions

   // Get value of normalVecId.
   template <int D> 
   inline int MaskGenFilmBase<D>::normalVecId() const
   {  return normalVecId_; }

   // Get value of interfaceThickness.
   template <int D> 
   inline double MaskGenFilmBase<D>::interfaceThickness() const
   {  return interfaceThickness_; }

   // Get value of excludedThickness.
   template <int D> 
   inline double MaskGenFilmBase<D>::excludedThickness() const
   {  return excludedThickness_; }

}
}
#endif