#ifndef PSPC_SOLVENT_H
#define PSPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>             // base class
#include <util/param/ParamComposite.h>     // base class
#include <pspc/solvers/Propagator.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf { 
namespace Pspc { 

   using namespace Util;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Solvent : public Species, public ParamComposite
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
      * Set association with Mesh and allocate concentration field array.
      *
      * \param mesh associated Mesh<D> object
      */
      void setMesh(Mesh<D> const & mesh);

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
      * \throw Exception if ensemble is open
      * \param phi desired chemical potential for this species
      */
      void setMu(double mu);

      /**
      * Set the monomer id for this solvent.
      *
      * \param monomerId integer id of monomer type (>=0)
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

      // Pointer to associated mesh
      Mesh<D> const *  meshPtr_;

   };
   
   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   template <int D>
   inline int Solvent<D>::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   template <int D>
   inline double Solvent<D>::size() const
   {  return size_; }
    
   /**
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
