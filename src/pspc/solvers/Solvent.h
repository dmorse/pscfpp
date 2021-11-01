#ifndef PSPC_SOLVENT_H
#define PSPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/field/RField.h>
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
      typedef RField<D> CField;

      /** 
      * Monomer chemical potential field type.
      */
      typedef RField<D> WField;

      /**
      * Constructor.
      */
      Solvent();
   
      /**
      * Constructor.
      */
      ~Solvent();
   
      /**
      * Set value of phi (volume fraction), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param phi  volume fraction for this species (input)
      */
      void setPhi(double phi);

      /**
      * Set value of mu (chemical potential), if ensemble is closed.
      *
      * \throw Exception if ensemble is closed
      * \param mu  chemical potential for this species (input)
      */
      void setMu(double mu);

      /**
      * Set association with Mesh and allocate concentration field array.
      *
      * \param mesh associated Mesh<D> object to describe spatial discretization.
      */
      void setMesh(Mesh<D> const & mesh);

      /**
      * Compute monomer concentration field and phi and/or mu.
      *
      * Upon return, concentration field, phi, mu, and q are all set.
      *
      * \param wField  chemical potential field of relevant monomer type (input)
      */
      void compute(WField const & wField );
   
      /**
      * Compute monomer concentration field and phi and/or mu.
      *
      * Takes array of wFields for all monomer types as an argument, but uses 
      * only the one with the correct monomer id for this solvent, by calling
      * compute(wFields[monomerId()]) internally. 
      *
      * \param wFields  array of monomer chemical potential fields (input)
      */ 
      void compute(DArray<WField> const & wFields);

      // Inherited functions

      using Pscf::SolventDescriptor::ensemble;

   protected:

      using Util::ParamComposite::setClassName;
      using Pscf::SolventDescriptor::size;
      using Pscf::SolventDescriptor::monomerId;

      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const
      {  return concentration_;  }
   
   private:

      CField concentration_;
  
      Mesh<D> const *  meshPtr_;
 
      using Pscf::Species::phi_;
      using Pscf::Species::mu_;

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
