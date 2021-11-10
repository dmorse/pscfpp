#ifndef PSPC_ITERATOR_H
#define PSPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(System<D>* system);

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

      /**
      * Is the unit cell adjusted to find a zero stress state?
      */
      bool isFlexible() const;

   protected:

      /// Flexible cell computation (true) or rigid (false), default value = false
      bool isFlexible_;

      /**
      * Get the parent system by reference.
      */
      System<D>& system();

   private:

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /**
      * Default constructor (private and not implemented to prohibit)
      */
      Iterator();

      /**
      * Copy constructor (private and not implemented to prohibit)
      */
      Iterator(System<D>& other);

   };

   template <int D>
   inline bool Iterator<D>::isFlexible() const
   { return isFlexible_; }

   template <int D>
   inline System<D>& Iterator<D>::system() 
   {  return *systemPtr_; }

} // namespace Pspc
} // namespace Pscf
#include "Iterator.tpp"
#endif
