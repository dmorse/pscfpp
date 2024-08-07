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
   * \ingroup Rpc_Iterator_Module
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

      using ParamComposite::setClassName;
      using MaskGenFilmBase<D>::normalVecId;
      using MaskGenFilmBase<D>::interfaceThickness;
      using MaskGenFilmBase<D>::excludedThickness;

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
      void setFlexibleParams();

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
      * Get the lattice parameters for this system.
      */
      FSArray<double, 6> systemLatticeParameters() const;

      /**
      * Get one of the lattice vectors for this system.
      * 
      * \param id  index of the desired lattice parameter
      */
      RealVec<D> systemLatticeVector(int id) const;

      using MaskGenFilmBase<D>::parametersCurrent_;

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
   {  return *sysPtr_; }

   // Get parent System by const reference.
   template <int D>
   System<D> const & MaskGenFilm<D>::system() const
   {  return *sysPtr_; }

   // Get space group name for this system.
   template <int D>
   inline std::string MaskGenFilm<D>::systemSpaceGroup() const
   {  return system().groupName(); }

   // Get lattice parameters for this system.
   template <int D>
   inline FSArray<double, 6> MaskGenFilm<D>::systemLatticeParameters() const
   {  return system().domain().unitCell().parameters(); }

   // Get one of the lattice vectors for this system.
   template <int D>
   inline RealVec<D> MaskGenFilm<D>::systemLatticeVector(int id) const
   {  return system().domain().unitCell().rBasis(id); }

}
}
#endif
