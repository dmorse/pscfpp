#ifndef PSPG_SOLVENT_H
#define PSPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventDescriptor.h>   // base class
#include <pspg/solvers/Propagator.h>       // typedefs
#include <pspg/field/RDField.h>

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
      * Compute monomer concentration field and phi and/or mu.
      *
      * Pure virtual function: Must be implemented by subclasses.
      * Upon return, concentration field, phi and mu are all set.
      *
      * \param wField monomer chemical potential field.
      */
      virtual void compute(WField const & wField );

      /**
      * Get the monomer concentration field for this solvent.
      */
      CField const & concField() const;
  
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
      CField concField_;
 
      /// Pointer to associated mesh
      Mesh<D> const *  meshPtr_;

   };

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D>
   inline const typename Solvent<D>::CField& Solvent<D>::concField() const
   {  return concField_;  }

}
} 
#endif 
