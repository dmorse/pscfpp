#ifndef RPG_MASK_GEN_FILM_H
#define RPG_MASK_GEN_FILM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/system/System.h>
#include <prdc/environment/FilmFieldGenMaskBase.h>  // Base class

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Field Generator for thin-film masks.
   *
   * The parent FilmFieldGenMaskBase class template defines all traits of a 
   * FilmFieldGenMask that do not require access to the System. This subclass
   * defines all methods that need System access.
   * 
   * If a MixAndMatchEnv contains a FilmFieldGenMask, then the system will
   * contain two parallel hard surfaces ("walls"), confining the
   * polymers/solvents to a "thin film" region of the unit cell.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class FilmFieldGenMask : public FilmFieldGenMaskBase<D>
   {

   public:

      /**
      * Default constructor
      */
      FilmFieldGenMask();
      
      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      FilmFieldGenMask(System<D>& sys);

      /**
      * Destructor
      */
      ~FilmFieldGenMask();

      /**
      * Get contribution to the stress from this mask
      * 
      * The mask defined by this class changes in a non-affine manner 
      * upon changing the lattice parameter corresponding to normalVecId.
      * Thus, if this lattice parameter is allowed to be flexible, the 
      * "stress" used to optimize the parameter must contain an additional 
      * contribution arising from the mask. This method evaluates this
      * contribution and returns its value. 
      * 
      * \param paramId  index of the lattice parameter being varied
      */
      double stress(int paramId) const;

      /**
      * Modify stress value in direction normal to the film.
      * 
      * The "stress" calculated by the System is used to minimize
      * fHelmholtz with respect to a given lattice parameter. In a thin 
      * film it is useful to instead minimize the excess free energy 
      * per unit area, (fHelmholtz - fRef) * Delta, where fRef is a 
      * reference free energy and Delta is the film thickness. The 
      * information needed to perform such a modification is contained 
      * within this object. This method performs this modification. The
      * stress will not be modified for lattice parameters that are 
      * parallel to the film.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      double modifyStress(int paramId, double stress) const;

      using FilmFieldGenMaskBase<D>::normalVecId;
      using FilmFieldGenMaskBase<D>::interfaceThickness;
      using FilmFieldGenMaskBase<D>::excludedThickness;
      using FilmFieldGenMaskBase<D>::hasFBulk;

   protected:

      /**
      * Compute the field and store where the System can access.
      */
      void compute();

      /**
      * Modifies iterator().flexibleParams_ to be compatible with the mask.
      */
      void setFlexibleParams() const;

      /**
      * Get the System associated with this object by reference.
      */
      System<D>& system();
      /**
      * Get the System associated with this object by const reference.
      */
      System<D> const & system() const;

      /**
      * Get the space group name for this system.
      */
      std::string systemSpaceGroup() const;

      /**
      * Get one of the lattice vectors for this system.
      * 
      * \param id  index of the desired lattice vector
      */
      RealVec<D> systemLatticeVector(int id) const;

      using FilmFieldGenMaskBase<D>::modifyFlexibleParams;
      using FilmFieldGenMaskBase<D>::normalVecCurrent_;
      using FilmFieldGenMaskBase<D>::fBulk_;
      using ParamComposite::setClassName;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Inline member functions

   // Get parent System by non-const reference.
   template <int D>
   System<D>& FilmFieldGenMask<D>::system() 
   {  return *sysPtr_; }

   // Get parent System by const reference.
   template <int D>
   System<D> const & FilmFieldGenMask<D>::system() const
   {  return *sysPtr_; }

   // Get space group name for this system.
   template <int D>
   inline std::string FilmFieldGenMask<D>::systemSpaceGroup() const
   {  return system().domain().groupName(); }

   // Get one of the lattice vectors for this system.
   template <int D>
   inline RealVec<D> FilmFieldGenMask<D>::systemLatticeVector(int id) const
   {  return system().domain().unitCell().rBasis(id); }

   // Explicit instantiation declarations
   extern template class FilmFieldGenMask<1>;
   extern template class FilmFieldGenMask<2>;
   extern template class FilmFieldGenMask<3>;

} // namespace Rpg
} // namespace Pscf

#endif
