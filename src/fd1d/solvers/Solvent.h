#ifndef FD1D_SOLVENT_H
#define FD1D_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>             // base class
#include <util/param/ParamComposite.h>     // base class
#include <fd1d/solvers/Propagator.h>
#include <fd1d/domain/Domain.h>

namespace Pscf { 
namespace Fd1d { 

   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ingroup Fd1d_Solver_Module
   */
   class Solvent : public Species, public ParamComposite
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
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

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
 
      /// \name Setters (set member data)
      //@{

      /**
      * Set value of phi (volume fraction), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param phi desired volume fraction for this species
      */
      void setPhi(double phi);

      /**
      * Set value of mu (chemical potential), if ensemble is closed.
      *
      * \throw Exception if ensemble is closed
      * \param mu  desired chemical potential for this species
      */
      void setMu(double mu);

      /**
      * Set the monomer id for this solvent.
      *
      * \param monomerId  integer id of monomer type, in [0,nMonomer-1]
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the size or volume of this solvent species.
      *
      * The ``size" is steric volume / reference volume.
      *
      * \param size volume of solvent
      */ 
      void setSize(double size);

      //@}
      /// \name Accessors (getters)
      //@{
 
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the size (number of monomers) in this solvent.
      */
      double size() const;

      /**
      * Return associated domain by reference.
      */
      Domain const & domain() const;

      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& cField() const;
  
      //@}

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

   private:

      // Concentration field for this solvent
      CField cField_;
  
      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume. 
      double size_;

      // Pointer to associated domain
      Domain const *  domainPtr_;

   };
   
   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   inline int Solvent::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   inline double Solvent::size() const
   {  return size_; }
    
   /**
   * Get monomer concentration field for this solvent.
   */
   inline const typename Solvent::CField& Solvent::cField() const
   {  return cField_;  }

   // Inline member functions

   /// Get Domain by reference.
   inline Domain const & Solvent::domain() const
   {   
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

}
}
#endif 
