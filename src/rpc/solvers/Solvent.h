#ifndef RPC_SOLVENT_H
#define RPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>   // base class
#include <prdc/cpu/RField.h>            // member

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Rpc { 

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref user_param_solvent_sec "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Solvent : public SolventSpecies
   {

   public:

      /**
      * Constructor.
      */
      Solvent();
   
      /**
      * Destructor.
      */
      ~Solvent();
  
      /**
      * Create an association with the mesh.
      *
      * \param mesh associated Mesh<D> object
      */
      void associate(Mesh<D> const & mesh);

      /**
      * Allocate memory for concentrationf field.
      */
      void allocate();

      /**
      * Compute monomer concentration field, q and phi and/or mu.
      *
      * Computes monomer concentration field cField, partition function
      * q, and either the solvent volume fraction phi or solvent chemical 
      * potential mu, depending on ensemble. The function takes the
      * chemical potential field wField for the relevant monomer type as 
      * its only input argument.
      *
      * The optional parameter phiTot is only relevant to problems such 
      * as thin films in which the material is excluded from part of the 
      * unit cell by imposing an inhomogeneous constraint on the sum of 
      * monomer concentrations (i.e., a "mask"). 
      *
      * \param wField  monomer chemical potential field of relevant type.
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(RField<D> const & wField, double phiTot = 1.0);

      /**
      * Get the monomer concentration field for this solvent.
      */
      RField<D> const & cField() const;
 
      // Inherited accessor functions 
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;
      using Pscf::SolventSpecies::monomerId;
      using Pscf::SolventSpecies::size;

   protected:

      // Inherited protected functions
      using Pscf::Species::setQ;
      
   private:

      /// Concentration field for this solvent.
      RField<D> cField_;
 
      /// Pointer to associated mesh.
      Mesh<D> const *  meshPtr_;

   };
   
   // Inline member function

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D>
   inline RField<D> const & Solvent<D>::cField() const
   {  return cField_;  }
  
   #ifndef RPC_SOLVENT_TPP
   // Supress implicit instantiation
   extern template class Solvent<1>;
   extern template class Solvent<2>;
   extern template class Solvent<3>;
   #endif

}
}
#endif 
