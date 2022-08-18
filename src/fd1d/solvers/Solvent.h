#ifndef FD1D_SOLVENT_H
#define FD1D_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventDescriptor.h>   // base class
//#include <pscf/chem/Species.h>             // base class
//#include <util/param/ParamComposite.h>     // base class
#include <fd1d/solvers/Propagator.h>
#include <fd1d/domain/Domain.h>

namespace Pscf { 
namespace Fd1d { 

   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref fd1d_Solvent_page "Parameter File Format"
   * \ingroup Fd1d_Solver_Module
   */
   //class Solvent : public Species, public ParamComposite
   class Solvent : public SolventDescriptor
   {

   public:

      /**
      * Monomer concentration field type.
      */
      typedef typename Propagator::CField CField;

      /** 
      * Monomer chemical potential field type.
      */
      typedef typename Propagator::WField WField;

      /**
      * Constructor.
      */
      Solvent();
   
      /**
      * Constructor.
      */
      ~Solvent();
  
      /**
      * Set association with Domain and allocate concentration field array.
      *
      * \param domain associated Domain object
      */
      void setDiscretization(Domain const & domain);

      /**
      * Compute monomer concentration field, q and phi and/or mu.
      *
      * Upon return, cField, phi, mu, and q are all set.
      *
      * \param wField monomer chemical potential field of relevant type.
      */
      void compute(WField const & wField );

      /// \name Accessors (getters)
      ///@{

      /**
      * Return associated domain by reference.
      */
      Domain const & domain() const;

      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& cField() const;
  
      ///@}

      // Inherited accessor functions 
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;

   protected:

      // Inherited data members
      using Pscf::Species::phi_;
      using Pscf::Species::mu_;
      using Pscf::Species::q_;
      using Pscf::Species::ensemble_;
      using Pscf::SolventDescriptor::monomerId_;
      using Pscf::SolventDescriptor::size_;

   private:

      // Concentration field for this solvent
      CField cField_;
 
      // Pointer to associated domain
      Domain const *  domainPtr_;

   };
   
   // Inline member functions

   /// Get monomer concentration field for this solvent.
   inline const typename Solvent::CField& Solvent::cField() const
   {  return cField_;  }

   /// Get associated Domain by reference.
   inline Domain const & Solvent::domain() const
   {   
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

}
}
#endif 
