#ifndef RPC_MASK_GEN_FILM_H
#define RPC_MASK_GEN_FILM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/System.h>
#include <prdc/iterator/MaskGenFilmBase.h>  // Base class

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Mask generator for a thin film geometry.
   *
   * The parent MaskGenFilmBase class template defines all traits of a 
   * MaskGenFilm that do not require access to the System. This subclass
   * defines all methods that need System access.
   * 
   * If the user chooses a MaskGenFilm as their mask generator, then the 
   * system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell.
   *
   * \ingroup Rpc_Scft_Iterator_Module
   */
   template <int D>
   class MaskGenFilm : public MaskGenFilmBase<D>
   {

   public:

      /**
      * Default constructor
      */
      MaskGenFilm();
      
      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      MaskGenFilm(System<D>& sys);

      /**
      * Destructor
      */
      ~MaskGenFilm();

      /**
      * Check whether the field has been generated
      */
      bool isGenerated() const;

      /**
      * Get contribution to the stress from this mask
      * 
      * The mask defined by this class changes in a non-affine manner 
      * upon changing the lattice parameter corresponding to normalVecId.
      * Thus, if this lattice parameter is allowed to be flexible, the 
      * "stress" used to optimize the parameter must contain additional 
      * terms arising from the mask. This method evaluates these terms
      * and returns their value. 
      * 
      * \param paramId  index of the lattice parameter being varied
      */
      double stressTerm(int paramId) const;

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
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      double modifyStress(int paramId, double stress) const;

      using MaskGenFilmBase<D>::normalVecId;
      using MaskGenFilmBase<D>::interfaceThickness;
      using MaskGenFilmBase<D>::excludedThickness;
      using MaskGenFilmBase<D>::hasFBulk;

   protected:

      /**
      * Allocate container necessary to generate and store field
      */ 
      void allocate();

      /**
      * Generate the field and store where the Iterator can access.
      */
      void generate();

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

      using MaskGenFilmBase<D>::modifyFlexibleParams;
      using MaskGenFilmBase<D>::normalVecCurrent_;
      using MaskGenFilmBase<D>::fBulk_;
      using ParamComposite::setClassName;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Inline member functions

   // Check whether the field has been generated
   template <int D>
   inline bool MaskGenFilm<D>::isGenerated() const
   {  return system().mask().hasData(); }

   // Get parent System by non-const reference.
   template <int D>
   System<D>& MaskGenFilm<D>::system() 
   {
      UTIL_CHECK(sysPtr_);  
      return *sysPtr_; 
   }

   // Get parent System by const reference.
   template <int D>
   System<D> const & MaskGenFilm<D>::system() const
   {  
      UTIL_CHECK(sysPtr_);  
      return *sysPtr_; 
   }

   // Get space group name for this system.
   template <int D>
   inline std::string MaskGenFilm<D>::systemSpaceGroup() const
   {  return system().domain().groupName(); }

   // Get one of the lattice vectors for this system.
   template <int D>
   inline RealVec<D> MaskGenFilm<D>::systemLatticeVector(int id) const
   {  return system().domain().unitCell().rBasis(id); }

   #ifndef RPC_MASK_GEN_FILM_TPP
   extern template class MaskGenFilm<1>;
   extern template class MaskGenFilm<2>;
   extern template class MaskGenFilm<3>;
   #endif

} // namespace Rpc
} // namespace Pscf

#endif
