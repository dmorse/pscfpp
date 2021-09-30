#ifndef PSPC_SOLVENT_H
#define PSPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/field/RField.h>
#include <pspc/solvers/Propagator.h>
#include <pscf/chem/Species.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/chem/SolventDescriptor.h>
#include <util/param/ParamComposite.h>

namespace Pscf { 
namespace Pspc { 

   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Solvent : public SolventDescriptor, public ParamComposite
   {
   public:

      /**
      * Monomer concentration field type.
      */
      typedef Propagator<D>::CField CField;

      /** 
      * Monomer chemical potential field type.
      */
      typedef Propagator<D>::WField WField;

      /**
      * Constructor.
      */
      Solvent();
   
      /**
      * Constructor.
      */
      ~Solvent();
   
      /**
      * Set association with Mesh and allocate concentration field array.
      *
      * \param mesh associated Mesh<D> object
      */
      void setMesh(Mesh<D> const & mesh);

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
      * \throw Exception if ensemble is open
      * \param phi desired chemical potential for this species
      */
      void setMu(double mu);

      /**
      * Compute monomer concentration field, q and phi and/or mu.
      *
      * Upon return, cField, phi, mu, and q are all set.
      *
      * \param wField monomer chemical potential field of relevant monomer type.
      */
      void solve(WField const & wField );
  
      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const
      {  return concentration_;  }
   
      // Inherited accessor functions 
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;
      using Pscf::SolventDescriptor::monomerId;
      using Pscf::SolventDescriptor::size;

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
      CField concentration_;
  
      // Pointer to associated mesh
      Mesh<D> const *  meshPtr_;
 
   };

   #ifndef PSPC_SOLVENT_TPP
   // Supress implicit instantiation
   extern template class Solvent<1>;
   extern template class Solvent<2>;
   extern template class Solvent<3>;
   #endif

}
} 
#endif 
