#ifndef PSPC_SOLVENT_H
#define PSPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventDescriptor.h>   // base class
#include <pspc/solvers/Propagator.h>       // typedefs

namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Pspc { 

   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref pspc_Solvent_page "Parameter File Format"
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   //class Solvent : public Species, public ParamComposite
   class Solvent : public SolventDescriptor
   {

   public:

      /**
      * Monomer concentration field type.
      */
      typedef typename Propagator<D>::CField CField;

      /** 
      * Monomer chemical potential field type.
      */
      typedef typename Propagator<D>::WField WField;

      /**
      * Constructor.
      */
      Solvent();
   
      /**
      * Destructor.
      */
      ~Solvent();
  
      /**
      * Set association with Mesh and allocate concentration field array.
      *
      * \param mesh associated Mesh<D> object
      */
      void setDiscretization(Mesh<D> const & mesh);

      /**
      * Compute monomer concentration field, q and phi and/or mu.
      *
      * Computes monomer concentration field cField, partition function
      * q, and either the solvent volume fraction phi or solvent chemical 
      * potential mu, depending on ensemble. The function takes the
      * chemical potential field wField for the relevant monomer type as 
      * its only input argument.
      *
      * \param wField  monomer chemical potential field of relevant type.
      */
      void compute(WField const & wField );

      /**
      * Get the monomer concentration field for this solvent.
      */
      CField const & cField() const;
 
      // Inherited accessor functions 
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;
      using Pscf::SolventDescriptor::monomerId;
      using Pscf::SolventDescriptor::size;

   protected:

      // Inherited protected data members
      using Pscf::Species::phi_;
      using Pscf::Species::mu_;
      using Pscf::Species::q_;
      using Pscf::Species::ensemble_;
      using Pscf::SolventDescriptor::monomerId_;
      using Pscf::SolventDescriptor::size_;

   private:

      /// Concentration field for this solvent
      CField cField_;
 
      /// Pointer to associated mesh
      Mesh<D> const *  meshPtr_;

   };
   
   // Inline member function

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D>
   inline const typename Solvent<D>::CField& Solvent<D>::cField() const
   {  return cField_;  }
  
   #ifndef PSPC_SOLVENT_TPP
   // Supress implicit instantiation
   extern template class Solvent<1>;
   extern template class Solvent<2>;
   extern template class Solvent<3>;
   #endif

}
}
#endif 
