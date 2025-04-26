#ifndef RPG_EXT_GEN_FILM_H
#define RPG_EXT_GEN_FILM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/System.h>
#include <prdc/iterator/ExtGenFilmBase.h>  // Base class

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * External field generator for a thin film geometry.
   *
   * The parent ExtGenFilmBase class template defines all traits of a 
   * ExtGenFilm that do not require access to the System. This subclass
   * defines all methods that need System access.
   * 
   * If the user chooses an ExtGenFilm object to generate external fields,
   * the external fields will have the same shape as the mask, with a 
   * magnitude defined by a Flory--Huggins-like chi parameter. This class
   * is specific to thin-film systems because it also allows for a 
   * different chi parameter to be defined on the top boundary than on
   * the bottom, through user input arrays chi_bottom and chi_top. See 
   * \ref scft_thin_films_page for more information. 
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class ExtGenFilm : public ExtGenFilmBase<D>
   {

   public:

      /**
      * Default constructor
      */
      ExtGenFilm();
      
      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      ExtGenFilm(System<D>& sys);

      /**
      * Destructor
      */
      ~ExtGenFilm();

      /**
      * Check whether the fields have been generated
      */
      bool isGenerated() const;

      /**
      * Get contribution to the stress from the external fields
      * 
      * The external fields defined by this class change in a non-affine 
      * manner upon changing the lattice parameter corresponding to 
      * normalVecId. Thus, if this lattice parameter is allowed to be 
      * flexible, the "stress" used to optimize the parameter must 
      * contain an additional term arising from the external fields. This 
      * method evaluates this term and returns its value. 
      * 
      * \param paramId  index of the lattice parameter being varied
      */
      double stressTerm(int paramId) const;

      using ExtGenFilmBase<D>::isAthermal;
      using ExtGenFilmBase<D>::chiBottom;
      using ExtGenFilmBase<D>::chiTop;
      using ExtGenFilmBase<D>::normalVecId;
      using ExtGenFilmBase<D>::interfaceThickness;
      using ExtGenFilmBase<D>::excludedThickness;

   protected:

      /**
      * Allocate container necessary to generate and store fields
      */ 
      void allocate();

      /**
      * Generate the fields and store where the Iterator can access.
      */
      void generate();

      /**
      * Get the System associated with this object by reference.
      */
      System<D> & system();
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

      /**
      * Get the number of monomer species for this system.
      */
      int systemNMonomer() const;

      /**
      * Use the mask to determine and store the value of normalVecId
      */
      void maskNormalVecId();

      /**
      * Use the mask to determine and store the value of interfaceThickness
      */
      void maskInterfaceThickness();

      using ExtGenFilmBase<D>::normalVecCurrent_;
      using ExtGenFilmBase<D>::chiBottomCurrent_;
      using ExtGenFilmBase<D>::chiTopCurrent_;
      using ParamComposite::setClassName;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

      /// interfaceThickness of the mask, obtained via maskInterfaceThickness
      double interfaceThickness_;

   };

   // Inline member functions

   // Check whether the field has been generated
   template <int D>
   inline bool ExtGenFilm<D>::isGenerated() const
   {  return system().h().hasData(); }

   // Get parent System by non-const reference.
   template <int D>
   System<D>& ExtGenFilm<D>::system() 
   {  return *sysPtr_; }

   // Get parent System by const reference.
   template <int D>
   System<D> const & ExtGenFilm<D>::system() const
   {  return *sysPtr_; }

   // Get space group name for this system.
   template <int D>
   inline std::string ExtGenFilm<D>::systemSpaceGroup() const
   {  return system().domain().groupName(); }

   // Get one of the lattice vectors for this system.
   template <int D>
   inline RealVec<D> ExtGenFilm<D>::systemLatticeVector(int id) const
   {  return system().domain().unitCell().rBasis(id); }

   // Get the number of monomer species for this system.
   template <int D>
   inline int ExtGenFilm<D>::systemNMonomer() const
   {  return system().mixture().nMonomer(); }

   #ifndef RPG_EXT_GEN_FILM_TPP
   extern template class ExtGenFilm<1>;
   extern template class ExtGenFilm<2>;
   extern template class ExtGenFilm<3>;
   #endif

} // namespace Rpg
} // namespace Pscf

#endif
