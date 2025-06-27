#ifndef RPC_EXT_GEN_FILM_H
#define RPC_EXT_GEN_FILM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/System.h>
#include <prdc/environment/FilmFieldGenExtBase.h>  // Base class

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Field Generator for external fields in thin-film systems.
   *
   * The parent FilmFieldGenExtBase class template defines all traits of a 
   * FilmFieldGenExt that do not require access to the System. This subclass
   * defines all methods that need System access.
   * 
   * If the user chooses an FilmFieldGenExt object to generate external 
   * fields, the external fields will have the same shape as the mask, 
   * with a magnitude defined by a Flory--Huggins-like chi parameter. This 
   * class is specific to thin-film systems because it also allows for a 
   * different chi parameter to be defined on the top boundary than on
   * the bottom, through user input arrays chi_bottom and chi_top. See 
   * \ref scft_thin_films_page for more information. 
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class FilmFieldGenExt : public FilmFieldGenExtBase<D>
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
      FilmFieldGenExt(System<D>& sys);

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

      using FilmFieldGenExtBase<D>::normalVecCurrent_;
      using FilmFieldGenExtBase<D>::chiBottomCurrent_;
      using FilmFieldGenExtBase<D>::chiTopCurrent_;
      using ParamComposite::setClassName;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

      /// Mask interfaceThickness, obtained via maskInterfaceThickness
      double interfaceThickness_;

   };

   // Inline member functions

   // Get parent System by non-const reference.
   template <int D>
   System<D>& FilmFieldGenExt<D>::system() 
   {
      UTIL_CHECK(sysPtr_);  
      return *sysPtr_; 
   }

   // Get parent System by const reference.
   template <int D>
   System<D> const & FilmFieldGenExt<D>::system() const
   {  
      UTIL_CHECK(sysPtr_);  
      return *sysPtr_; 
   }

   // Get space group name for this system.
   template <int D>
   inline std::string FilmFieldGenExt<D>::systemSpaceGroup() const
   {  return system().domain().groupName(); }

   // Get one of the lattice vectors for this system.
   template <int D>
   inline RealVec<D> FilmFieldGenExt<D>::systemLatticeVector(int id) const
   {  return system().domain().unitCell().rBasis(id); }

   // Get the number of monomer species for this system.
   template <int D>
   inline int FilmFieldGenExt<D>::systemNMonomer() const
   {  return system().mixture().nMonomer(); }
   
   #ifndef RPC_EXT_GEN_FILM_TPP
   extern template class FilmFieldGenExt<1>;
   extern template class FilmFieldGenExt<2>;
   extern template class FilmFieldGenExt<3>;
   #endif

} // namespace Rpc
} // namespace Pscf

#endif
