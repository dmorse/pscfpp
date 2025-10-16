#ifndef RPG_SOLVENT_H
#define RPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>   // base class
#include <prdc/cuda/RField.h>           // member

// Forward declaration
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf { 
namespace Rpg { 

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref user_param_solvent_sec "Manual Page"
   * \ingroup Rpg_Solver_Module
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
      * Associate this object with a mesh.
      * 
      * Must be called before allocate().
      *
      * \param mesh  Mesh<D> object - spatial discretization mesh
      */
      void associate(Mesh<D> const & mesh);

      /**
      * Allocate internal data containers. 
      * 
      * associate() must have been called first.
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
 
   private:

      /// Concentration field for this solvent
      RField<D> cField_;
 
      /// Pointer to associated mesh
      Mesh<D> const *  meshPtr_;

   };

   /*
   * Associate this object with a mesh.
   */
   template <int D>
   inline void Solvent<D>::associate(Mesh<D> const & mesh)
   {
      UTIL_CHECK(mesh.size() > 1);
      meshPtr_ = &mesh;
   }

   /*
   * Allocate internal data containers. 
   */
   template <int D>
   inline void Solvent<D>::allocate()
   {  
      UTIL_CHECK(meshPtr_);
      cField_.allocate(meshPtr_->dimensions()); 
   }

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D>
   inline RField<D> const & Solvent<D>::cField() const
   {  return cField_;  }

   // Explicit instantiation declarations
   extern template class Solvent<1>;
   extern template class Solvent<2>;
   extern template class Solvent<3>;

}
} 
#endif 
