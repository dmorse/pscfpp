#ifndef RPC_EXT_GEN_FILM_H
#define RPC_EXT_GEN_FILM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/System.h>
#include <prdc/iterator/ExtGenFilmBase.h>  // Base class

namespace Pscf {
namespace Rpc {

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
   * \ingroup Rpc_Iterator_Module
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

      using ExtGenFilmBase<D>::isAthermal;
      using ExtGenFilmBase<D>::chiBottom;
      using ExtGenFilmBase<D>::chiTop;
      using ParamComposite::setClassName;

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
      * Get the lattice parameters for this system.
      */
      FSArray<double, 6> systemLatticeParameters() const;

      /**
      * Get the number of monomer species for this system.
      */
      int systemNMonomer() const;

      /**
      * Use the mask to determine the value of normalVecId
      */
      void maskNormalVecId();

      using ExtGenFilmBase<D>::parametersCurrent_;
      using ExtGenFilmBase<D>::chiBottomCurrent_;
      using ExtGenFilmBase<D>::chiTopCurrent_;
      using ExtGenFilmBase<D>::normalVecId_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

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
   {  return system().groupName(); }

   // Get lattice parameters for this system.
   template <int D>
   inline FSArray<double, 6> ExtGenFilm<D>::systemLatticeParameters() const
   {  return system().domain().unitCell().parameters(); }

   // Get the number of monomer species for this system.
   template <int D>
   inline int ExtGenFilm<D>::systemNMonomer() const
   {  return system().mixture().nMonomer(); }

}
}
#endif
