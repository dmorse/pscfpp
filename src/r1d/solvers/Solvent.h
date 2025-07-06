#ifndef R1D_SOLVENT_H
#define R1D_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>   // base class
#include <r1d/solvers/Propagator.h>     // typedefs

namespace Pscf { 
namespace R1d { 

   class Domain;
   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref user_param_solvent_sec "Parameter File Format"
   * \ingroup R1d_Solver_Module
   */
   class Solvent : public SolventSpecies
   {

   public:

      // Public typename alias

      /**
      * Field type.
      */
      using FieldT = typename Propagator::FieldT;

      // Public member functions

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
      void compute(FieldT const & wField );

      /// \name Accessors (getters)
      ///@{

      /**
      * Return associated domain by reference.
      */
      Domain const & domain() const;

      /**
      * Get monomer concentration field for this solvent.
      */
      const FieldT& cField() const;
  
      ///@}

   private:

      // Concentration field for this solvent
      FieldT cField_;
 
      // Pointer to associated domain
      Domain const *  domainPtr_;

   };
   
   // Inline member functions

   /// Get monomer concentration field for this solvent.
   inline const typename Solvent::FieldT& Solvent::cField() const
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
