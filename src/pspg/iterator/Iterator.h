#ifndef PSPG_ITERATOR_H
#define PSPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspg_Iterator_Module
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
      Iterator(System<D>* system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve() = 0;

      System<D>* systemPtr_;
   };

} // namespace Pspg
} // namespace Pscf
#include "Iterator.tpp"
#endif
