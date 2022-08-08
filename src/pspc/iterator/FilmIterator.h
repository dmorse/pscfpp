#ifndef PSPC_FILM_ITERATOR_H
#define PSPC_FILM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIteratorBase.h"

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * Descriptor for a FilmIterator object. The parent FilmIteratorBase 
   * class template defines all traits of a FilmIterator that do not 
   * depend on D, the dimension of the system. This FilmIterator class 
   * defines only the partial specializations of FilmIterator for 1D, 
   * 2D, and 3D.
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

   template <int D>
   class System;

   template <int D, typename IteratorType>
   class FilmIterator : public FilmIteratorBase<D,IteratorType>
   {};

   // Partial Specializations

   // 1D
   template <typename IteratorType>
   class FilmIterator<1,IteratorType> : public FilmIteratorBase<1,IteratorType>
   {
   public:

      /**
      * Constructor.
      */
      FilmIterator(System<1>& system);

      /**
      * Output the indices of each flexible lattice parameter, based on
      * normalVec and unitCell definitions in param file. Assumes that
      * isFlexible == true, and gives a warning if none of the parameters
      * are actually able to be varied given the thin film constraints
      */
      FSArray<int, 6> flexibleParams() const;

      /**
      * Check that user-defined lattice basis vectors (stored in the
      * Domain<D> object associated with this FilmIterator class)
      * are compatible with the thin film constraint
      */
      void checkLatticeVectors() const;

      using FilmIteratorBase<1,IteratorType>::setClassName;
      using FilmIteratorBase<1,IteratorType>::normalVec;

   protected:

      using FilmIteratorBase<1,IteratorType>::system;
      using FilmIteratorBase<1,IteratorType>::iterator;
      
   };

   // 2D
   template <typename IteratorType>
   class FilmIterator<2,IteratorType> : public FilmIteratorBase<2,IteratorType>
   {
   public:

      /**
      * Constructor.
      */
      FilmIterator(System<2>& system);

      /**
      * Output the indices of each flexible lattice parameter, based on
      * normalVec and unitCell definitions in param file. Assumes that
      * isFlexible == true, and gives a warning if none of the parameters
      * are actually able to be varied given the thin film constraints
      */
      FSArray<int, 6> flexibleParams() const;

      /**
      * Check that user-defined lattice basis vectors (stored in the
      * Domain<D> object associated with this FilmIterator class)
      * are compatible with the thin film constraint
      */
      void checkLatticeVectors() const;

      using FilmIteratorBase<2,IteratorType>::setClassName;
      using FilmIteratorBase<2,IteratorType>::normalVec;

   protected:

      using FilmIteratorBase<2,IteratorType>::system;
      using FilmIteratorBase<2,IteratorType>::iterator;
      
   };

   // 3D
   template <typename IteratorType>
   class FilmIterator<3,IteratorType> : public FilmIteratorBase<3,IteratorType>
   {
   public:

      /**
      * Constructor.
      */
      FilmIterator(System<3>& system);

      /**
      * Output the indices of each flexible lattice parameter, based on
      * normalVec and unitCell definitions in param file. Assumes that
      * isFlexible == true, and gives a warning if none of the parameters
      * are actually able to be varied given the thin film constraints
      */
      FSArray<int, 6> flexibleParams() const;

      /**
      * Check that user-defined lattice basis vectors (stored in the
      * Domain<D> object associated with this FilmIterator class)
      * are compatible with the thin film constraint
      */
      void checkLatticeVectors() const;
      
      using FilmIteratorBase<3,IteratorType>::setClassName;
      using FilmIteratorBase<3,IteratorType>::normalVec;

   protected:

      using FilmIteratorBase<3,IteratorType>::system;
      using FilmIteratorBase<3,IteratorType>::iterator;
      
   };

} // namespace Pspc
} // namespace Pscf
#endif
