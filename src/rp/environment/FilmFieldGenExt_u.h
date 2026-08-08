#ifndef RPG_EXT_GEN_FILM_U_H
#define RPG_EXT_GEN_FILM_U_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/environment/FilmFieldGenExtBase.h>  // base class template
#include <pscf/backends/CUT.h>                      // base class argument
#include <rpg/system/System.h>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class FilmFieldGenExt;
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Field Generator for external fields in thin-film systems.
   *
   * The parent FilmFieldGenExtBase class template defines all traits of
   * a FilmFieldGenExt that do not require access to the System. This 
   * subclass defines all methods that need System access.
   * 
   * If the user chooses a FilmFieldGenExt object to generate external 
   * fields, each external fields will have the same dependence on the
   * normal coordinate  as the mask, with a prefactor defined by a 
   * Flory--Huggins-like chi parameter.  This class is specific to 
   * thin-film systems because it also allows for a different chi 
   * parameter to be defined on the top boundary than on the bottom, 
   * through user input arrays chi_bottom and chi_top. See 
   * \ref scft_thin_films_page for more information. 
   *
   * \ingroup Rp_Field_Module
   */
   template <int D>
   class FilmFieldGenExt<D,CUT> 
     : public FilmFieldGenExtBase<D>
   {

   public:

      /**
      * Default constructor
      */
      FilmFieldGenExt();
      
      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      FilmFieldGenExt(System<D,CUT>& sys);

      /**
      * Destructor
      */
      ~FilmFieldGenExt();

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
      double stress(int paramId) const;

      using FilmFieldGenExtBase<D>::isAthermal;
      using FilmFieldGenExtBase<D>::chiBottom;
      using FilmFieldGenExtBase<D>::chiTop;
      using FilmFieldGenExtBase<D>::normalVecId;
      using FilmFieldGenExtBase<D>::interfaceThickness;
      using FilmFieldGenExtBase<D>::excludedThickness;

   protected:

      /**
      * Compute the fields and store where the System can access.
      */
      void compute();

      /**
      * Get the System associated with this object by reference.
      */
      System<D,CUT> & system();

      /**
      * Get the System associated with this object by const reference.
      */
      System<D,CUT> const & system() const;

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

      using FilmFieldGenExtBase<D>::normalVecCurrent_;
      using FilmFieldGenExtBase<D>::chiBottomCurrent_;
      using FilmFieldGenExtBase<D>::chiTopCurrent_;
      using ParamComposite::setClassName;

   private:

      /// Pointer to the associated system object.
      System<D,CUT>* sysPtr_;

      /// interfaceThickness of the mask, obtained via maskInterfaceThickness
      double interfaceThickness_;

   };

   // Inline member functions

   // Get parent System by non-const reference.
   template <int D> inline
   System<D,CUT>& 
   FilmFieldGenExt<D,CUT>::system() 
   {  return *sysPtr_; }

   // Get parent System by const reference.
   template <int D> inline
   System<D,CUT> const & 
   FilmFieldGenExt<D,CUT>::system() const
   {  return *sysPtr_; }

   // Get space group name for this system.
   template <int D> inline 
   std::string 
   FilmFieldGenExt<D,CUT>::systemSpaceGroup() const
   {  return system().domain().groupName(); }

   // Get one of the lattice vectors for this system.
   template <int D> inline 
   RealVec<D> 
   FilmFieldGenExt<D,CUT>::systemLatticeVector(int id) const
   {  return system().domain().unitCell().rBasis(id); }

   // Get the number of monomer species for this system.
   template <int D> inline 
   int FilmFieldGenExt<D,CUT>::systemNMonomer() const
   {  return system().mixture().nMonomer(); }

   // Explicit instantiation declarations
   extern template class FilmFieldGenExt<1,CUT>;
   extern template class FilmFieldGenExt<2,CUT>;
   extern template class FilmFieldGenExt<3,CUT>;

} // namespace Rp
} // namespace Pscf
#endif
