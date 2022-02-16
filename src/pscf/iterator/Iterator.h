#ifndef PSCF_ITERATOR_H
#define PSCF_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include "IteratorMediator.h"
#include <util/global.h>                  

namespace Pscf {

   template <typename T>
   class IteratorMediator;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pscf_Iterator_Module
   */

   template <typename T>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(IteratorMediator<T>& iterMed);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Setup and allocate memory before iteration.
      */
      virtual void setup() = 0;

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve() = 0;

   protected:

      /**
      * Get the iterator mediator by reference.
      */
      IteratorMediator<T>& iterMed();

   private:

      /// Pointer to iterator mediator object
      IteratorMediator<T>* iterMed_;

      /**
      * Default constructor (private and not implemented to prohibit)
      */
      Iterator();

      /**
      * Copy constructor (private and not implemented to prohibit)
      */
      Iterator(Iterator& other);

   };

   template <typename T>
   inline IteratorMediator<T>& Iterator<T>::iterMed() 
   {  return *iterMed_; }

   template <typename T>
   Iterator<T>::Iterator(IteratorMediator<T>& iterMed)
    : iterMed_(&iterMed)
   {  setClassName("Iterator"); }

   template <typename T>
   Iterator<T>::~Iterator()
   {}

} // namespace Pscf
#endif
