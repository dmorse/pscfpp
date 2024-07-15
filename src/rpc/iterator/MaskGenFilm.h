#ifndef RPC_MASK_GEN_FILM_H
#define RPC_MASK_GEN_FILM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/iterator/MaskGenFilmBase.h>  // Base class
#include <rpc/System.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Mask generator for a thin film geometry (empty base template).
   *
   * The parent MaskGenFilmBase class template defines all traits of a 
   * MaskGenFilm that do not depend on D, the dimension of space. This 
   * MaskGenFilm class template is an empty template that is replaced
   * by partial specializations for D=1, 2 and 3.
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
      * \param itr  Iterator parent object
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
      void generateMask();

      /**
      * Modifies iterator().flexibleParams_ to be compatible with the mask.
      */
      void setFlexibleParams();

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
      std::string getSpaceGroup() const;

      /**
      * Get the lattice parameters for this system.
      */
      FSArray<double, 6> getLatticeParameters() const;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

      using MaskGenFilmBase<D>::parameters_;

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
   inline std::string MaskGenFilm<D>::getSpaceGroup() const
   {  return system().groupName(); }

   // Get lattice parameters for this system.
   template <int D>
   inline FSArray<double, 6> MaskGenFilm<D>::getLatticeParameters() const
   {  return system().domain().unitCell().parameters(); }

}
}
#endif
