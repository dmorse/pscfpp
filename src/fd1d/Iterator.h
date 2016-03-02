#ifndef FD1D_ITERATOR_H
#define FD1D_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Fd1d
{

   class Mixture;
   class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Iterator();

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Create association with the parent System.
      * 
      * \param system parent System object.
      */
      virtual void setSystem(System& system);

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve() = 0;

   protected:

      System& system();
      
      Mixture& mixture();
      
   private:

      // Pointer to parent System object.
      System* systemPtr_;

      // Pointer to associated Mixture object.
      Mixture* mixturePtr_;

   };

   inline
   System& Iterator::system()
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   inline
   Mixture& Iterator::mixture()
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
