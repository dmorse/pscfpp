#ifndef PSPC_FILM_ITERATOR_H
#define PSPC_FILM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include "FilmIteratorBase.h"

namespace Pscf {
namespace Rpc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Iterator for a thin film (empty base template).
   *
   * The parent FilmIteratorBase class template defines all traits of a 
   * FilmIterator that do not depend on D, the dimension of space. This 
   * FilmIterator class template is an empty template that is replaced
   * by partial specializations for D=1, 2 and 3.
   * 
   * If the user chooses a FilmIterator as their iterator, then the 
   * system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell.
   * This only affects the iterator, not the rest of SCFT, so we isolate
   * the imposition of the thin film constraint to this subclass of 
   * iterator. This is essentially a wrapper for any other type of 
   * iterator (e.g., AmIterator), that adds the additional functionality 
   * required to impose the thin film constraint properly.
   *
   * FilmIterator is generalized to be compatible with any iterator within
   * it, as long as the iterator can impose 1) a mask that confines the 
   * polymers/solvents to a certain region of space, and 2) an external 
   * field. A FilmIterator object in the param file is created by appending
   * "Film" to the end of the name of the iterator that is stored inside of
   * FilmIterator (e.g., "AmIteratorFilm{" will create a FilmIterator with
   * an AmIterator object inside of it). 
   *
   * \ingroup Rpc_Iterator_Module
   */
   template <int D, typename IteratorType>
   class FilmIterator : public FilmIteratorBase<D,IteratorType>
   {};

   // Partial Specializations

   /**
   * FilmIterator specialization for 1D problems.
   */
   template <typename IteratorType>
   class FilmIterator<1, IteratorType> 
      : public FilmIteratorBase<1, IteratorType>
   {
   public:

      /**
      * Constructor.
      */
      FilmIterator(System<1>& system);

      /**
      * Modifies flexibleParams_ to be compatible with thin film constraint.
      *
      * Construct an array indicating whether each lattice parameter is 
      * flexible, based on normalVecId and unitCell definitions in param 
      * file as well as the optional user input flexibleParams. Store this
      * array in flexibleParams_ member of this object, as well the 
      * flexibleParams_ member of the iterator within this object.
      * Uses the flexibleParams_ member of the iterator within this object
      * as a starting point.
      * 
      * In 1D, a thin film can not have flexible lattice parameters, so 
      * this will always set flexibleParams_ to an array of zeroes.
      */
      void setFlexibleParams();

      /**
      * Check compatibility of lattice with thin film constraint.
      *
      * Check that the user-defined lattice basis vectors in the
      * associated Domain<D> object are compatible with the thin 
      * film constraint
      */
      void checkLatticeVectors() const;

      using FilmIteratorBase<1,IteratorType>::setClassName;
      using FilmIteratorBase<1,IteratorType>::normalVecId;

   protected:

      using FilmIteratorBase<1,IteratorType>::system;
      using FilmIteratorBase<1,IteratorType>::iterator;
      using Iterator<1>::setFlexibleParams;
      
   };

   /**
   * FilmIterator specialization for 2D problems (confined to strip).
   */
   template <typename IteratorType>
   class FilmIterator<2,IteratorType> 
      : public FilmIteratorBase<2,IteratorType>
   {
   public:

      /**
      * Constructor.
      */
      FilmIterator(System<2>& system);

      /**
      * Modifies flexibleParams_ to be compatible with thin film constraint.
      * 
      * Construct an array indicating whether each lattice parameter is
      * flexible, based on normalVecId and unitCell definitions in param 
      * file as well as the optional user input flexibleParams. Store this
      * array in flexibleParams_ member of this object, as well the 
      * flexibleParams_ member of the iterator within this object.
      * Uses the flexibleParams_ member of the iterator within this object
      * as a starting point.
      */
      void setFlexibleParams();

      /**
      * Check compatibility of lattice with thin film constraint.
      *
      * Check that the user-defined lattice basis vectors in the
      * associated Domain<D> object are compatible with the thin 
      * film constraint
      */
      void checkLatticeVectors() const;

      using FilmIteratorBase<2,IteratorType>::setClassName;
      using FilmIteratorBase<2,IteratorType>::normalVecId;
      using Iterator<2>::nFlexibleParams;

   protected:

      using FilmIteratorBase<2,IteratorType>::system;
      using FilmIteratorBase<2,IteratorType>::iterator;
      using Iterator<2>::setFlexibleParams;
      
   };

   /**
   * FilmIterator specialization for 3D problems (confined to slit).
   */
   template <typename IteratorType>
   class FilmIterator<3,IteratorType> 
      : public FilmIteratorBase<3,IteratorType>
   {
   public:

      /**
      * Constructor.
      */
      FilmIterator(System<3>& system);

      /**
      * Modifies flexibleParams_ to be compatible with thin film constraint.
      * 
      * Construct an array indicating whether each lattice parameter is
      * flexible, based on normalVecId and unitCell definitions in param 
      * file as well as the optional user input flexibleParams. Store this
      * array in flexibleParams_ member of this object, as well the 
      * flexibleParams_ member of the iterator within this object.
      * Uses the flexibleParams_ member of the iterator within this object
      * as a starting point.
      */
      void setFlexibleParams();

      /**
      * Check compatibility of lattice with thin film constraint.
      *
      * Check that the user-defined lattice basis vectors in the
      * associated Domain<D> object are compatible with the thin 
      * film constraint
      */
      void checkLatticeVectors() const;
      
      using FilmIteratorBase<3,IteratorType>::setClassName;
      using FilmIteratorBase<3,IteratorType>::normalVecId;
      using Iterator<3>::nFlexibleParams;

   protected:

      using FilmIteratorBase<3,IteratorType>::system;
      using FilmIteratorBase<3,IteratorType>::iterator;
      using Iterator<3>::setFlexibleParams;
      
   };

   #ifndef PSPC_FILM_ITERATOR_TPP
   // Suppress implicit instantiation
   extern template class FilmIterator<1, AmIteratorBasis<1> >;
   extern template class FilmIterator<2, AmIteratorBasis<2> >;
   extern template class FilmIterator<3, AmIteratorBasis<3> >;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif
