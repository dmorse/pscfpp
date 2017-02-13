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
   class Domain;
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
      * \param isContinuation true iff part of sweep, and not first step.
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve(bool isContinuation = false) = 0;

   protected:

      System & system();

      System const & system() const;
      
      Mixture& mixture();

      Mixture const & mixture() const;
      
      Domain& domain();

      Domain const & domain() const;
      
   private:

      // Pointer to parent System object.
      System* systemPtr_;

      // Pointer to associated Mixture object.
      Mixture* mixturePtr_;

      // Pointer to associated Domain object.
      Domain* domainPtr_;

   };

   inline
   System& Iterator::system()
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   inline
   System const & Iterator::system() const
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

   inline
   Mixture const & Iterator::mixture() const
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

   inline
   Domain& Iterator::domain()
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   inline
   Domain const & Iterator::domain() const
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
