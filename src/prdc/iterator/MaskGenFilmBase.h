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
   * scft_thin_films_page for more information. 
   * 
   * \ingroup Prdc_Iterator_Module
   */
   template <int D>
   class MaskGenFilmBase : public FieldGenerator
   {

   public:

      /**
      * Constructor.
      */
      MaskGenFilmBase();

      /**
      * Destructor.
      */
      ~MaskGenFilmBase();

      /**
      * Read parameter file block and initialize.
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
      * Get value of fBulk.
      */
      double fBulk() const;

      /**
      * Check whether a value of fBulk was provided.
      */
      bool hasFBulk() const;

      /**
      * Check whether the field has been generated.
      */
      bool isGenerated() const = 0;

      /**
      * Get contribution to the stress from this mask
      * 
      * The mask defined by this class changes in a non-affine manner 
      * upon changing the lattice parameter corresponding to normalVecId.
      * Thus, if this lattice parameter is allowed to be flexible, the 
      * "stress" used to optimize the parameter must contain additional 
      * terms arising from the mask. This method evaluates these terms
      * and returns their value. Access to various System properties is
      * required, so the method must be implemented by subclasses.
      * 
      * \param paramId  index of the lattice parameter being varied
      */
      virtual double stressTerm(int paramId) const = 0;

      /**
      * Modify stress value in direction normal to the film.
      * 
      * The "stress" calculated by the Mixture object is used to minimize
      * fHelmholtz with respect to a given lattice parameter. In a thin 
      * film it is useful to instead minimize the excess free energy 
      * per unit area, (fHelmholtz - fRef) * Delta, where fRef is a 
      * reference free energy and Delta is the film thickness. The 
      * information needed to perform such a modification is contained 
      * within this object. This method performs this modification. The
      * stress will not be modified for lattice parameters that are 
      * parallel to the film.
      * 
      * This method requires access to various attributes of the System,
      * so it is left as a pure virtual method here and should be 
      * implemented by subclasses.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      virtual double modifyStress(int paramId, double stress) const = 0;

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
      * constant. Subclasses should define this method to set the 
      * flexibility of these lattice parameters appropriately so that 
      * the film thickness is held constant, as are all but one of the 
      * angles between basis vectors (the angle in the plane of the film
      * may vary). The other lattice parameters are allowed to be 
      * flexible if the user specified that they are flexible in the 
      * Iterator parameters. 
      * 
      * Note that the lattice parameter that defines the film thickness
      * may be flexible if the optional input parameter fBulk is provided.
      * This parameter is necessary to compute the stress in the direction
      * normal to the film.
      */
      virtual void setFlexibleParams() = 0;

      /**
      * Get the space group name for this system.
      */
      virtual std::string systemSpaceGroup() const = 0;

      /**
      * Get one of the lattice vectors for this system.
      * 
      * \param id  index of the desired lattice vector
      */
      virtual RealVec<D> systemLatticeVector(int id) const = 0;

      /**
      * The lattice vector normal to the film used to generate these fields.
      * 
      * This vector is set to be equal to the system's lattice vector with
      * index normalVecId_ each time the external fields are generated. The 
      * system's lattice vectors may then change, and this normalVecCurrent_
      * vector is used to detect whether they have changed. This is used to 
      * decide whether a new set of external fields needs to be generated.
      */
      RealVec<D> normalVecCurrent_;

      /// Reference free energy used to calculate stress normal to the film
      double fBulk_;

      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Lattice basis vector that is normal to the walls
      int normalVecId_;

      /// Interface thickness
      double interfaceThickness_;

      /// Excluded (wall) thickness
      double excludedThickness_;

      /// Does this object have a value of fBulk from the parameter file?
      bool hasFBulk_;

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

   // Get value of fBulk.
   template <int D> 
   inline double MaskGenFilmBase<D>::fBulk() const
   {  return fBulk_; }

   // Check whether a value of fBulk was provided.
   template <int D> 
   inline bool MaskGenFilmBase<D>::hasFBulk() const
   {  return hasFBulk_; }

}
}
#endif
