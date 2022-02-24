#ifndef PSPC_ITERATOR_H
#define PSPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <pspc/System.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspc {

   template <int D>
   class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspc_Iterator_Module
   */

   template <int D>
   class Iterator : virtual public ParamComposite
   {

   public:

      /**
      * Constructor.
      * 
      * \param system system object by reference
      */
      Iterator(System<D>& system);

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
      * Get the system by reference.
      */
      System<D>& system();

      /// Pointer to system object
      System<D>* sys_;

   private:

      /**
      * Default constructor (private and not implemented to prohibit)
      */
      Iterator();

      /**
      * Copy constructor (private and not implemented to prohibit)
      */
      Iterator(Iterator& other);

   };

   template <int D>
   inline System<D>& Iterator<D>::system() 
   {  return *sys_; }

   template <int D>
   Iterator<D>::Iterator(System<D>& system)
    : sys_(&system)
   {  setClassName("Iterator"); }

   template <int D>
   Iterator<D>::~Iterator()
   {}

} // namespace Pspc
} // namespace Pscf
#endif
