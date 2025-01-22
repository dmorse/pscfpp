#ifndef RPG_SOLVENT_H
#define RPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventDescriptor.h>   // base class
#include <rpg/solvers/Propagator.h>       // typedefs
#include <prdc/cuda/RField.h>

namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf { 
namespace Rpg { 

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ingroup Rpg_Solvers_Module
   */
   template <int D>
   class Solvent : public SolventDescriptor
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
      * Compute monomer concentration field and phi and/or mu.
      *
      * Upon return, concentration field, phi and mu are all set.
      *
      * \param wField monomer chemical potential field
      */
      void compute(RField<D> const & wField );

      /**
      * Get the monomer concentration field for this solvent.
      */
      RField<D> const & concField() const;
  
      // Inherited public accessor functions 
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
      RField<D> concField_;
 
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
      concField_.allocate(meshPtr_->dimensions()); 
   }

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D>
   inline RField<D> const & Solvent<D>::concField() const
   {  return concField_;  }

}
} 
#endif 
