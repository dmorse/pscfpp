#ifndef PSPC_FILM_ITERATOR_BASE_H
#define PSPC_FILM_ITERATOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "pspc/iterator/Iterator.h"        // base class
#include "util/containers/FSArray.h"       // container
#include <string>
#include <iostream>

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;
   
   using namespace Util;

   /**
   * Descriptor for a FilmIterator object. This base class template defines
   * all traits of a FilmIterator that do not depend on D, the dimension of
   * the system. The subclasses of FilmIteratorBase are the partial
   * specializations of FilmIterator for 1D, 2D, and 3D.
   * 
   * If the user chooses a FilmIterator as their iterator, then the 
   * system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell.
   * This only affects the iterator, not the rest of SCFT, so we isolate
   * the imposition of the thin film constraint to this subclass of 
   * iterator. This is essentially a wrapper for any other type of iterator
   * (e.g., AmIterator), that adds the additional functionality required
   * to impose the thin film constraint properly.
   * 
   * FilmIterator is generalized to be compatible with any iterator within
   * it, as long as the iterator can impose 1) a mask that confines the 
   * polymers/solvents to a certain region of space, and 2) an external 
   * field. A FilmIterator object in the param file is created by appending
   * "Film" to the end of the name of the iterator that is stored inside of
   * FilmIterator (e.g., "AmIteratorFilm{" will create a FilmIterator with
   * an AmIterator object inside of it). 
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D, typename IteratorType>
   class FilmIteratorBase : public Iterator<D>
   {

   public:

      /**
      * Constructor.
      */
      FilmIteratorBase(System<D>& system);

      /**
      * Destructor.
      */
      ~FilmIteratorBase();

      /**
      * Read and initialize.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Iterate to a solution
      *
      * \param isContinuation true iff continuation within a sweep
      */
      int solve(bool isContinuation = false);

      /**
      * Return const reference to the real iterator within this FilmIterator
      */
      IteratorType const & iterator() const;

      /**
      * Modifies flexibleParams_ to be compatible with thin film constraint.
      * 
      * Modifies the flexibleParams_ array to be compatible with the thin
      * film constraint. Stores resulting array in flexibleParams_ member
      * of this object, as well as the flexibleParams_ member of the
      * iterator within this object.
      * 
      * This function varies depending on D, the dimensionality of the
      * system. Therefore, it is implemented in the partial class 
      * specializations in FilmIterator, rather than in this base class.
      */
      virtual void setFlexibleParams() = 0;

      /**
      * Check that lattice vectors are compatible with thin film constraint.
      * 
      * Check that user-defined lattice basis vectors (stored in the
      * Domain<D> object associated with this FilmIterator class)
      * are compatible with the thin film constraint. All lattice basis
      * vectors must be either parallel or perpendicular to the walls.
      * 
      * This function varies depending on D, the dimensionality of the
      * system. Therefore, it is implemented in the partial class 
      * specializations in FilmIterator, rather than in this base class.
      */
      virtual void checkLatticeVectors() const = 0;

      /**
      * Generates mask and external field for the walls and stores in System.
      * 
      * Generates the field representation of the walls, based on the values
      * of wallThickness and interfaceThickness that were input by the user.
      * Then, stores this wall field in system().mask() to be used
      * as a mask during iteration, and also passes the corresponding
      * external potential fields into system().h() if isAthermal() = false.
      */
      void generateWallFields();

      /**
      * Updates the mask and external fields for the walls if needed.
      * 
      * Checks whether the lattice parameters have been updated since the
      * last call of generateWallFields(), and if the parameters have
      * changed then calls generateWallFields() again to update them. 
      * 
      * Also updates the external fields if the wall/polymer chi parameters
      * have been updated since external fields were last generated.
      */
      void updateWallFields();

      /**
      * Check that space group is compatible with the thin film constraint.
      */
      void checkSpaceGroup() const;

      /**
      * Are the walls chemically identical?
      * 
      * This is the case when chiBottom is equal to chiTop.
      */
      bool isSymmetric() const;

      /**
      * Are the walls athermal?
      * 
      * This is only true if all values in chiBottom and chiTop are zero.
      */
      bool isAthermal() const;

      /**
      * Set the value of chi between species s and the bottom wall.
      * 
      * \param s  species index, 0 <= id < nVertex
      * \param chi  value of chi(s,w)
      */
      void setChiBottom(int s, double chi);

      /**
      * Set the value of chi between species s and the top wall.
      * 
      * \param s  species index, 0 <= id < nVertex
      * \param chi  value of chi(s,w)
      */
      void setChiTop(int s, double chi);

      /**
      * Get value of normalVecId
      */
      int normalVecId() const;

      /**
      * Get value of interfaceThickness
      */
      double interfaceThickness() const;

      /**
      * Get value of wallThickness
      */
      double wallThickness() const;

      /**
      * Get const chiBottom matrix by reference
      */
      DArray<double> const & chiTop() const;

      /**
      * Get const chiTop array by reference
      */
      DArray<double> const & chiBottom() const;

      /**
      * Get the chi parameter between the bottom wall and species s
      * 
      * \param s  species index, 0 <= id < nVertex
      */
      double chiBottom(int s) const;

      /**
      * Get the chi parameter between the top wall and species s
      * 
      * \param s  species index, 0 <= id < nVertex
      */
      double chiTop(int s) const;

      using Iterator<D>::isFlexible;

   protected:

      /**
      * Initialize just before entry to iterative loop.
      * 
      * Allocate required memory, perform necessary checks to ensure user 
      * input is compatible with a film constraint, and create the mask / 
      * external fields that will be used to represent the walls during 
      * iteration.
      */
      void setup();
      
      /**
      * Return reference to the real iterator within this FilmIterator
      */
      IteratorType& iterator();

      /**
      * Generate external fields for the walls.
      * 
      * Generate external fields only, and pass them into System for use
      * during calculation. This is called by generateWallFields().
      */
      void generateExternalFields();
       
      void outputTimers(std::ostream& out){};
      void clearTimers(){};
      using Iterator<D>::system;
      using Iterator<D>::setClassName;
      using Iterator<D>::isFlexible_;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::readDArray;
      using ParamComposite::setParent;
      using ParamComposite::addComponent;

   private:

      /// The actual iterator that does all the work
      IteratorType iterator_;

      /// Lattice parameters associated with the current maskBasis
      FSArray<double, 6> parameters_;

      /// Lattice basis vector that is normal to the walls
      int normalVecId_;

      /// Interface thickness
      double t_;

      /// Wall thickness
      double T_;

      /// chiBottom array
      DArray<double> chiBottom_;

      /// chiTop array
      DArray<double> chiTop_;

      /// Wall chiBottom array associated with the current system().h() field
      DArray<double> chiBottomCurrent_;

      /// Wall chiTop array associated with the current system().h() field
      DArray<double> chiTopCurrent_;

      /// Flag indicating whether the wall fields are currently ungenerated
      bool ungenerated_;
   };

   // Inline member functions

   // Return reference to iterator within this FilmIterator
   template <int D, typename IteratorType>
   inline IteratorType& FilmIteratorBase<D, IteratorType>::iterator()
   {  return iterator_; }

   // Return const reference to iterator within this FilmIterator
   template <int D, typename IteratorType>
   inline 
   IteratorType const & FilmIteratorBase<D, IteratorType>::iterator() const
   {  return iterator_; }

   // Set value of chi between species s and the bottom wall
   template <int D, typename IteratorType>
   inline 
   void FilmIteratorBase<D, IteratorType>::setChiBottom(int s, double chi)
   {  chiBottom_[s] = chi; }

   // Set value of chi between species s and the top wall
   template <int D, typename IteratorType>
   inline 
   void FilmIteratorBase<D, IteratorType>::setChiTop(int s, double chi)
   {  chiTop_[s] = chi; }

   // Get value of normalVecId
   template <int D, typename IteratorType>
   inline int FilmIteratorBase<D, IteratorType>::normalVecId() const
   {  return normalVecId_; }

   // Get value of interfaceThickness
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::interfaceThickness() 
   const
   {  return t_; }

   // Get value of wallThickness
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::wallThickness() const
   {  return T_; }

   // Get chiBottom array by const reference
   template <int D, typename IteratorType>
   inline 
   DArray<double> const & FilmIteratorBase<D, IteratorType>::chiBottom() 
   const
   {  return chiBottom_; }

   // Get chiTop array by const reference
   template <int D, typename IteratorType>
   inline 
   DArray<double> const & FilmIteratorBase<D, IteratorType>::chiTop() 
   const
   {  return chiTop_; }

   // Get the chi parameter between the bottom wall and species s
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::chiBottom(int s) 
   const
   {  return chiBottom_[s]; }

   // Get the chi parameter between the top wall and species s
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::chiTop(int s) 
   const
   {  return chiTop_[s]; }

} // namespace Pspc
} // namespace Pscf

#include "FilmIteratorBase.tpp"
#endif
