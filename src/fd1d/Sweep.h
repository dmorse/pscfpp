#ifndef FD1D_SWEEP_H
#define FD1D_SWEEP_H

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
   class Iterator;

   using namespace Util;

   /**
   * Base class for classes solve along a line in parameter space.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Sweep : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Sweep();

      /**
      * Destructor.
      */
      ~Sweep();

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

      System & system();

      System const & system() const;
      
      Mixture& mixture();

      Mixture const & mixture() const;
      
      Domain& domain();

      Domain const & domain() const;
      
      Iterator& iterator();

      Iterator const & iterator() const;
      
   private:

      // Pointer to parent System object.
      System* systemPtr_;

      // Pointer to associated Mixture object.
      Mixture* mixturePtr_;

      // Pointer to associated Domain object.
      Domain* domainPtr_;

      // Pointer to associated Iterator object.
      Iterator* iteratorPtr_;

   };

   inline
   System& Sweep::system()
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   inline
   System const & Sweep::system() const
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   inline
   Mixture& Sweep::mixture()
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

   inline
   Mixture const & Sweep::mixture() const
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

   inline
   Domain& Sweep::domain()
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   inline
   Domain const & Sweep::domain() const
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   inline
   Iterator& Sweep::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   inline
   Iterator const & Sweep::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
