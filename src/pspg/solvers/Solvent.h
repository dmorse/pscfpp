#ifndef PSPG_SOLVENT_H
#define PSPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventDescriptor.h>   // base class
#include <pspg/solvers/Propagator.h>       // typedefs
#include <pspg/field/RField.h>

namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf { 
namespace Pspg { 

   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ingroup Pspg_Solvers_Module
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
      * Set association with Mesh, allocate memory.
      *
      * \param mesh associated Mesh<D> object (input)
      */
      void setDiscretization(Mesh<D> const & mesh);
   
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
   * Get monomer concentration field for this solvent.
   */
   template <int D>
   inline RField<D> const & Solvent<D>::concField() const
   {  return concField_;  }

}
} 
#endif 
