#ifndef PRDC_FILM_MASK_FG_BASE_H
#define PRDC_FILM_MASK_FG_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>            // Function parameter
#include <pscf/environment/FieldGenerator.h>  // Base class
#include <pscf/math/RealVec.h>                // container
#include <iostream>
#include <string>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base class Field Generator for thin-film masks.
   * 
   * This is a base class for FilmFieldGenMask that defines all traits of a 
   * FilmFieldGenMask that do not require access to the System (System access is
   * needed, for example, to get the space group and set the mask field).
   * 
   * If the user chooses a FilmFieldGenMask object to construct the mask, then 
   * the system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell. The
   * shape of the mask is defined by three input parameters: normalVecId,
   * excludedThickness, and interfaceThickness. See \ref 
   * scft_thin_films_page for more information. 
   * 
   * \ingroup Prdc_Field_Module
   */
   template <int D>
   class FilmFieldGenMaskBase : public FieldGenerator
   {

   public:

      /**
      * Constructor.
      */
      FilmFieldGenMaskBase();

      /**
      * Destructor.
      */
      ~FilmFieldGenMaskBase();

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
      virtual double stress(int paramId) const = 0;

      /**
      * Modify stress value in direction normal to the film.
      * 
      * The "stress" calculated by the System object is used to minimize
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
      void checkLatticeVectors() const;

      /**
      * Allocate container necessary to generate and store field.
      */ 
      virtual void allocate() = 0;

      /**
      * Generate the field and store where the System can access.
      */
      virtual void generate() = 0;

      /**
      * Sets iterator's flexibleParams array to be compatible with the mask.
      * 
      * If the iterator allows for flexible lattice parameters, this 
      * method will access the "flexibleParams" array of the iterator
      * (a FSArray<bool,6> that indicates which parameters are flexible),
      * create a modified version using modifyFlexibleParams to be 
      * compatible with this mask, and then update the flexibleParams 
      * array owned by the iterator to match this modified array.
      * 
      * Because this requires access to the iterator, it must be 
      * implemented by subclasses.
      */
      virtual void setFlexibleParams() const = 0;

      /**
      * Modifies a flexibleParams array to be compatible with this mask.
      * 
      * An iterator that is compatible with this mask should either have 
      * a strictly rigid unit cell (lattice parameters are fixed), or 
      * should allow some lattice parameters to be flexible while others
      * are rigid. In the latter case, an FSArray<bool,6> will be used
      * to keep track of which parameters are flexible and which are not.
      *  
      * This "flexibleParams" array may not be compatible with the 
      * constraints of this mask, so this method will read the current
      * array, modify it for compatibility, and then return the corrected
      * version. In particular, this method will ensure that the lattice
      * parameter defining the film thickness is rigid, unless the 
      * optional input parameter fBulk is provided (since this parameter
      * is necessary to compute the stress in the direction normal to
      * the film). Further, all angles between basis vectors should be 
      * fixed except the angle in the plane of the film.
      * 
      * Subclasses should use this method in generate() to set the
      * flexibleParams array of the iterator.
      * 
      * \param current flexibleParams array 
      * \param cell unit cell (read only)
      */
      FSArray<bool,6> modifyFlexibleParams(FSArray<bool,6> current,
                                           UnitCell<D> const & cell) const;

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
      * index normalVecId_ each time the mask is generated. The system's
      * lattice vectors may then change, and this normalVecCurrent_ vector
      * is used to detect whether they have changed. This is used to decide
      * whether a new mask needs to be generated.
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
   inline int FilmFieldGenMaskBase<D>::normalVecId() const
   {  return normalVecId_; }

   // Get value of interfaceThickness.
   template <int D> 
   inline double FilmFieldGenMaskBase<D>::interfaceThickness() const
   {  return interfaceThickness_; }

   // Get value of excludedThickness.
   template <int D> 
   inline double FilmFieldGenMaskBase<D>::excludedThickness() const
   {  return excludedThickness_; }

   // Get value of fBulk.
   template <int D> 
   inline double FilmFieldGenMaskBase<D>::fBulk() const
   {  return fBulk_; }

   // Check whether a value of fBulk was provided.
   template <int D> 
   inline bool FilmFieldGenMaskBase<D>::hasFBulk() const
   {  return hasFBulk_; }

   #ifndef PRDC_MASK_GEN_FILM_BASE_TPP
   // Suppress implicit instantiation
   extern template class FilmFieldGenMaskBase<1>;
   extern template class FilmFieldGenMaskBase<2>;
   extern template class FilmFieldGenMaskBase<3>;
   #endif

}
}
#endif
