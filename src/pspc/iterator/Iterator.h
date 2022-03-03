#ifndef PSPC_ITERATOR_H
#define PSPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/DArray.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   typedef DArray<double> FieldCPU;

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
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(System<D>& system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Setup iterator.
      */
      virtual void setup() = 0;

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve() = 0;

      /**
      * Return true if unit cell is flexible, false if rigid.
      */
      bool isFlexible() 
      {  return isFlexible_;}

      /**
      * Return reference to parent system.
      */
      System<D>& system() 
      {  return *sysPtr_;}

   protected:

      /// Is the unit cell flexible during iteration?
      bool isFlexible_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;
      
   };

   template<int D>
   inline Iterator<D>::Iterator()
   {  setClassName("Iterator"); }

   template<int D>
   inline Iterator<D>::Iterator(System<D>& system)
   : sysPtr_(&system)
   {  setClassName("Iterator"); }

   template<int D>
   inline Iterator<D>::~Iterator()
   {}

} // namespace Pspc
} // namespace Pscf
#endif
